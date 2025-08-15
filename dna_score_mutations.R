# Required libraries
library(msa)
library(Biostrings)

# ----- Core scorer (DNA) -----
# Reference-anchored indexing + robust single-base deletion detection (±1 around mapped col).
calculate_mutation_score <- function(reference, target, changes_str) {
  if (is.null(changes_str) || changes_str == "" || is.na(changes_str)) {
    return(list(score = 0, detected_mutations = "None",
                wt_confirmed_positions = 0, status = "Unknown"))
  }

  # Parse mutations like A198G, A198- (single base deletion)
  change_pairs <- strsplit(changes_str, ",")[[1]]
  change_pairs <- trimws(change_pairs)

  change_data <- data.frame(
    Position     = numeric(0),
    ReferenceDNA = character(0),
    TargetDNA    = character(0),
    stringsAsFactors = FALSE
  )

  for (chg in change_pairs) {
    if (nchar(chg) < 3) next
    pos <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chg)))
    if (is.na(pos)) next
    ref_base <- toupper(substr(chg, 1, 1))
    tgt_base <- toupper(substr(chg, nchar(chg), nchar(chg)))  # "-" for deletion
    change_data <- rbind(change_data, data.frame(
      Position = pos, ReferenceDNA = ref_base, TargetDNA = tgt_base
    ))
  }

  tryCatch({
    # Normalize sequences: uppercase, strip whitespace; IMPORTANT: no pre-gapped reads
    reference <- toupper(gsub("[ \r\n\t]", "", reference))
    target    <- toupper(gsub("[- \r\n\t]", "", target))
    
    # Pad the TARGET (read) with N's to prevent terminal-gap artefacts in msa()
    pad_len <- 4000L
    pad_blk <- paste(rep("N", pad_len), collapse = "")
    target  <- paste0(pad_blk, target, pad_blk)

    # Align with ClustalW
    sequences <- DNAStringSet(c(Reference = reference, Target = target))
    aln <- msa(sequences, method = "ClustalW", verbose = FALSE)

    # Extract aligned strings (with gaps)
    aligned_ref    <- as.character(unmasked(aln)["Reference"])
    aligned_target <- as.character(unmasked(aln)["Target"])
    ref_chars    <- strsplit(aligned_ref, "")[[1]]
    target_chars <- strsplit(aligned_target, "")[[1]]

    # Map ungapped reference positions -> alignment columns
    ref_pos_to_col <- which(ref_chars != "-")
    max_ref_pos <- length(ref_pos_to_col)

    score <- 0
    detected_mutations <- character(0)
    wt_confirmed <- 0

    for (i in seq_len(nrow(change_data))) {
      pos        <- change_data$Position[i]
      ref_DNA    <- change_data$ReferenceDNA[i]
      target_DNA <- change_data$TargetDNA[i]  # "-" for deletion, or base for substitution

      if (is.na(pos) || pos < 1 || pos > max_ref_pos) next

      col_idx <- ref_pos_to_col[pos]  # alignment col for that reference position
      base_ref_at_pos  <- ref_chars[col_idx]
      base_read_at_pos <- target_chars[col_idx]

      # --- Deletion handling (e.g., A198-) ---
      if (target_DNA == "-") {
        # check column, column-1, column+1 to tolerate 1-col jitter
        called_del <- FALSE

        check_gap_at <- function(j) {
          if (j < 1 || j > length(target_chars)) return(FALSE)
          # ensure mutation spec matches the reference residue at pos
          if (ref_chars[col_idx] != ref_DNA) return(FALSE)
          # count deletion if read shows gap at j and ref at j is a real base (not '-')
          target_chars[j] == "-" && ref_chars[j] != "-"
        }

        if (check_gap_at(col_idx) ||
            check_gap_at(col_idx - 1) ||
            check_gap_at(col_idx + 1)) {
          score <- score + 1
          detected_mutations <- c(detected_mutations, paste0(ref_DNA, pos, "-"))
          called_del <- TRUE
        }

        # If no deletion called but explicit WT at mapped column, confirm WT
        if (!called_del && base_read_at_pos == ref_DNA) {
          wt_confirmed <- wt_confirmed + 1
        }
        next
      }

      # --- Substitutions & WT confirmation ---
      # Treat gap/unknown in READ as uninformative
      if (base_read_at_pos %in% c("-", "N")) next

      if (base_read_at_pos == target_DNA) {
        # Resistance substitution detected
        score <- score + 1
        detected_mutations <- c(detected_mutations, paste0(ref_DNA, pos, target_DNA))
      } else if (base_ref_at_pos == ref_DNA && base_read_at_pos == ref_DNA) {
        # Explicit WT confirmation
        wt_confirmed <- wt_confirmed + 1
      }
      # else: other mismatch -> ignore under simplified rule
    }

    detected_mutations_str <- if (length(detected_mutations)) paste(detected_mutations, collapse = ", ") else "None"

    # Status rule: Resistant if any; else Wildtype if any WT; else Unknown
    status <- if (score > 0) "Resistant" else if (wt_confirmed > 0) "Wildtype" else "Unknown"

    list(
      score = score,
      detected_mutations = detected_mutations_str,
      wt_confirmed_positions = wt_confirmed,
      status = status
    )

  }, error = function(e) {
    warning(paste("Error processing sequences:", e$message))
    list(score = NA, detected_mutations = "Error",
         wt_confirmed_positions = NA, status = "Error")
  })
}

# ----- Runner: reads TSV and writes annotated output -----
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

    results <- lapply(seq_len(nrow(input_data)), function(i) {
      calculate_mutation_score(
        reference   = input_data$reference[i],
        target      = input_data$target[i],
        changes_str = input_data$changes_str[i]
      )
    })

    input_data$MutationScore         <- sapply(results, `[[`, "score")
    input_data$DetectedMutations     <- sapply(results, `[[`, "detected_mutations")
    input_data$WTConfirmedPositions  <- sapply(results, `[[`, "wt_confirmed_positions")
    input_data$Status                <- sapply(results, `[[`, "status")

    write.table(input_data,
                "updated_table_with_scores_and_mutations.tsv",
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)

    input_data
  }, error = function(e) {
    stop(paste("Error processing file:", e$message))
  })
}

result <- process_mutation_data("input.tsv")
print(result)
