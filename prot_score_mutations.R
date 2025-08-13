# Required libraries
library(msa)
library(Biostrings)

# Function to calculate mutation score and detect mutations (protein)
# Uses simple leading-offset model; calls substitutions & single-AA deletions
calculate_mutation_score <- function(reference, target, changes_str) {
  if (is.null(changes_str) || changes_str == "" || is.na(changes_str)) {
    return(list(score = 0, detected_mutations = "None", wt_confirmed_positions = 0, status = "Unknown"))
  }

  # Parse "A51V", "Q88-" etc.
  change_pairs <- strsplit(changes_str, ",")[[1]]
  change_pairs <- trimws(change_pairs)

  change_data <- data.frame(
    Position   = numeric(0),
    ReferenceAA = character(0),
    TargetAA    = character(0),
    stringsAsFactors = FALSE
  )

  for (change in change_pairs) {
    if (nchar(change) < 3) next
    position <- as.numeric(gsub("[^0-9]", "", change))
    if (is.na(position)) next
    reference_aa <- toupper(substr(change, 1, 1))
    target_aa    <- toupper(substr(change, nchar(change), nchar(change)))  # may be "-" for deletion
    change_data <- rbind(change_data,
                         data.frame(Position = position,
                                    ReferenceAA = reference_aa,
                                    TargetAA = target_aa))
  }

  tryCatch({
    # Normalize sequences (uppercase, strip whitespace)
    reference <- toupper(gsub("[ \r\n\t]", "", reference))
    target    <- toupper(gsub("[ \r\n\t]", "", target))

    # Align
    sequences <- AAStringSet(c(Reference = reference, Target = target))
    alignment <- msa(sequences, method = "ClustalW", verbose = FALSE)

    # Extract aligned strings (with gaps)
    aligned_ref    <- as.character(unmasked(alignment)["Reference"])
    aligned_target <- as.character(unmasked(alignment)["Target"])

    ref_chars    <- strsplit(aligned_ref, "")[[1]]
    target_chars <- strsplit(aligned_target, "")[[1]]

    # Simple fix: compute leading offset from first non-gap in reference
    first_ref_col <- which(ref_chars != "-")[1]
    offset <- if (length(first_ref_col)) first_ref_col - 1 else 0

    score <- 0
    detected_mutations <- character(0)
    wt_confirmed <- 0

    for (i in 1:nrow(change_data)) {
      pos <- change_data$Position[i]
      ref_AA <- change_data$ReferenceAA[i]
      target_AA <- change_data$TargetAA[i]  # may be "-"

      aln_idx <- offset + pos
      if (aln_idx < 1 || aln_idx > length(target_chars)) next

      aa_in_read <- target_chars[aln_idx]
      aa_in_ref  <- ref_chars[aln_idx]

      # Deletion call: expected deletion (TargetAA == "-") and read shows a gap at this ref column
      if (target_AA == "-") {
        if (aa_in_read == "-") {
          score <- score + 1
          detected_mutations <- c(detected_mutations, paste0(ref_AA, pos, "-"))
        }
        next
      }

      # Treat read '-' or 'X' as uninformative (neither WT nor R)
      if (aa_in_read %in% c("-", "X")) next

      # Resistance substitution?
      if (aa_in_read == target_AA) {
        score <- score + 1
        detected_mutations <- c(detected_mutations, paste0(ref_AA, pos, target_AA))
      # Explicit WT confirmation (ensure the reference at this column is indeed the ref AA)
      } else if (aa_in_ref == ref_AA && aa_in_read == ref_AA) {
        wt_confirmed <- wt_confirmed + 1
      }
      # else: other mismatch -> ignore for simplified logic
    }

    detected_mutations_str <- if (length(detected_mutations) > 0) {
      paste(detected_mutations, collapse = ", ")
    } else {
      "None"
    }

    # Simplified status logic
    if (score > 0) {
      status <- "Resistant"
    } else if (wt_confirmed > 0) {
      status <- "Wildtype"
    } else {
      status <- "Unknown"
    }

    return(list(
      score = score,
      detected_mutations = detected_mutations_str,
      wt_confirmed_positions = wt_confirmed,
      status = status
    ))

  }, error = function(e) {
    warning(paste("Error processing sequences:", e$message))
    return(list(score = NA, detected_mutations = "Error", wt_confirmed_positions = NA, status = "Error"))
  })
}

# Read/process the input file and write output
process_mutation_data <- function(file_path) {
  tryCatch({
    input_data <- read.delim(file_path,
                             sep = "\t",
                             stringsAsFactors = FALSE,
                             quote = "",
                             check.names = FALSE)

    required_cols <- c("class", "gene", "read", "reference", "target", "changes_str")
    missing_cols <- required_cols[!required_cols %in% colnames(input_data)]
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    }

    results <- lapply(1:nrow(input_data), function(i) {
      calculate_mutation_score(
        input_data$reference[i],
        input_data$target[i],
        input_data$changes_str[i]
      )
    })

    input_data$MutationScore         <- sapply(results, function(x) x$score)
    input_data$DetectedMutations     <- sapply(results, function(x) x$detected_mutations)
    input_data$WTConfirmedPositions  <- sapply(results, function(x) x$wt_confirmed_positions)
    input_data$Status                <- sapply(results, function(x) x$status)

    write.table(input_data,
                "updated_table_with_scores_and_mutations.tsv",
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)

    return(input_data)

  }, error = function(e) {
    stop(paste("Error processing file:", e$message))
  })
}

# Example usage:
result <- process_mutation_data("input.tsv")
print(result)
