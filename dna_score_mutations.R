# Required libraries
library(Biostrings)
library(pwalign)


calculate_mutation_score <- function(reference, target, changes_str, debug = FALSE) {
  if (is.null(changes_str) || changes_str == "" || is.na(changes_str)) {
    return(list(score = 0, detected_mutations = "None",
                wt_confirmed_positions = 0, status = "Unknown"))
  }

  # Parse mutations like A15G, G43-, etc.
  change_pairs <- strsplit(changes_str, ",")[[1]]
  change_pairs <- trimws(change_pairs)
  change_data <- data.frame(
    Position = numeric(0),
    ReferenceDNA = character(0),
    TargetDNA = character(0),
    stringsAsFactors = FALSE
  )
  for (chg in change_pairs) {
    if (nchar(chg) < 3) next
    pos <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chg)))
    if (is.na(pos)) next
    ref_base <- toupper(substr(chg, 1, 1))
    tgt_base <- toupper(substr(chg, nchar(chg), nchar(chg)))  # "-" for deletion
    change_data <- rbind(change_data,
                         data.frame(Position = pos, ReferenceDNA = ref_base, TargetDNA = tgt_base))
  }

  tryCatch({
    # Normalize inputs
    reference <- toupper(gsub("[ \r\n\t]", "", reference))
    target    <- toupper(gsub("[- \r\n\t]", "", target))

    # Align: force full reference (pattern) to align; read (subject) can overhang
    submat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -5, baseOnly = TRUE)
    pa <-  pairwiseAlignment(pattern = DNAString(reference),
                            subject = DNAString(target),
                            type = "global-local",
                            substitutionMatrix = submat,
                            gapOpening = -10, gapExtension = -0.5)

    aligned_ref    <- as.character(alignedPattern(pa))
    aligned_target <- as.character(alignedSubject(pa))
    ref_chars    <- strsplit(aligned_ref, "")[[1]]
    target_chars <- strsplit(aligned_target, "")[[1]]

    # Build alignment-column -> reference-position map (NA where ref has a gap)
    ref_pos_at_col <- rep(NA_integer_, length(ref_chars))
    rp <- 0L
    for (k in seq_along(ref_chars)) {
      if (ref_chars[k] != "-") {
        rp <- rp + 1L
        ref_pos_at_col[k] <- rp
      }
    }
    max_ref_pos <- rp

    informative <- function(b) b %in% c("A","C","G","T")

    score <- 0L
    wt_confirmed <- 0L
    detected <- character(0)

    for (i in seq_len(nrow(change_data))) {
      pos        <- change_data$Position[i]
      ref_DNA    <- change_data$ReferenceDNA[i]
      target_DNA <- change_data$TargetDNA[i]

      if (is.na(pos) || pos < 1 || pos > max_ref_pos) next

      # find the unique alignment column that corresponds to reference position `pos`
      cols <- which(ref_pos_at_col == pos)
      if (length(cols) != 1L) {
        if (debug) message(sprintf("Reference pos %d not uniquely mapped (cols: %s)", pos, paste(cols, collapse=",")))
        next
      }
      j <- cols[1]

      ref_at_j  <- ref_chars[j]
      read_at_j <- target_chars[j]

      # sanity: mutation spec must match the reference base at this position
      if (ref_at_j != ref_DNA) {
        if (debug) message(sprintf("Spec/ref mismatch at pos %d: spec=%s ref=%s", pos, ref_DNA, ref_at_j))
        next
      }

      if (target_DNA == "-") {
        # Deletion: read must have a gap exactly at the reference-mapped column
        if (read_at_j == "-" && ref_at_j != "-") {
          score <- score + 1L
          detected <- c(detected, paste0(ref_DNA, pos, "-"))
          if (debug) message(sprintf("DEL %s at pos %d (col %d)", paste0(ref_DNA, pos, "-"), pos, j))
        } else if (informative(read_at_j) && read_at_j == ref_DNA) {
          wt_confirmed <- wt_confirmed + 1L
          if (debug) message(sprintf("WT confirm (no deletion) at pos %d", pos))
        } else if (debug) {
          message(sprintf("No DEL at pos %d: read(col)=%s", pos, read_at_j))
        }
      } else {
        # Substitution: read base must equal target_DNA at the mapped column
        if (informative(read_at_j) && read_at_j == target_DNA) {
          score <- score + 1L
          detected <- c(detected, paste0(ref_DNA, pos, target_DNA))
          if (debug) message(sprintf("SUB %s at pos %d (col %d)", paste0(ref_DNA, pos, target_DNA), pos, j))
        } else if (informative(read_at_j) && read_at_j == ref_DNA) {
          wt_confirmed <- wt_confirmed + 1L
          if (debug) message(sprintf("WT confirm at pos %d", pos))
        } else if (debug) {
          message(sprintf("No SUB at pos %d: read(col)=%s target=%s", pos, read_at_j, target_DNA))
        }
      }
    }

    status <- if (score > 0) "Resistant" else if (wt_confirmed > 0) "Wildtype" else "Unknown"
    list(score = score,
         detected_mutations = if (length(detected)) paste(detected, collapse=", ") else "None",
         wt_confirmed_positions = wt_confirmed,
         status = status)

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
warnings()
print(result)
