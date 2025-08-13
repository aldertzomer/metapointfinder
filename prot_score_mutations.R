# Required libraries
library(msa)
library(Biostrings)

# ----- Core scorer (protein) -----
# Reference-anchored indexing + robust single-AA deletion detection (±1 around the mapped col).
calculate_mutation_score <- function(reference, target, changes_str) {
  if (is.null(changes_str) || changes_str == "" || is.na(changes_str)) {
    return(list(score = 0, detected_mutations = "None", wt_confirmed_positions = 0, status = "Unknown"))
  }

  # Parse mutations like A51V, Q42- (single AA del)
  change_pairs <- strsplit(changes_str, ",")[[1]]
  change_pairs <- trimws(change_pairs)

  change_data <- data.frame(
    Position    = numeric(0),
    ReferenceAA = character(0),
    TargetAA    = character(0),
    stringsAsFactors = FALSE
  )

  for (chg in change_pairs) {
    if (nchar(chg) < 3) next
    pos <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chg)))
    if (is.na(pos)) next
    ref_aa <- toupper(substr(chg, 1, 1))
    tgt_aa <- toupper(substr(chg, nchar(chg), nchar(chg)))  # "-" for deletion
    change_data <- rbind(change_data, data.frame(
      Position = pos, ReferenceAA = ref_aa, TargetAA = tgt_aa
    ))
  }

  tryCatch({
    # Normalize sequences (uppercase, strip whitespace). IMPORTANT: strip any literal gaps from raw target
    reference <- toupper(gsub("[ \r\n\t]", "", reference))
    target    <- toupper(gsub("[- \r\n\t]", "", target))

    # Align. Using ClustalW method as that does not delete flanks
    sequences <- AAStringSet(c(Reference = reference, Target = target))
    aln <- msa(sequences, method = "ClustalW", verbose = FALSE)

    # Extract aligned strings (with gaps preserved)
    aligned_ref    <- as.character(unmasked(aln)["Reference"])
    aligned_target <- as.character(unmasked(aln)["Target"])
    ref_chars    <- strsplit(aligned_ref, "")[[1]]
    target_chars <- strsplit(aligned_target, "")[[1]]

    # Reference-anchored: map ungapped reference positions -> alignment columns
    ref_pos_to_col <- which(ref_chars != "-")
    max_ref_pos <- length(ref_pos_to_col)

    score <- 0
    detected_mutations <- character(0)
    wt_confirmed <- 0

    for (i in seq_len(nrow(change_data))) {
      pos       <- change_data$Position[i]
      ref_AA    <- change_data$ReferenceAA[i]
      target_AA <- change_data$TargetAA[i]  # "-" for deletion, or a residue for substitution

      if (is.na(pos) || pos < 1 || pos > max_ref_pos) next

      col_idx <- ref_pos_to_col[pos]  # alignment column corresponding to reference position `pos`
      aa_ref_at_pos  <- ref_chars[col_idx]
      aa_read_at_pos <- target_chars[col_idx]

      # --- Deletions (e.g., Q42-) ---
      if (target_AA == "-") {
        # tolerate ClustalW left/right gap jitter by checking col_idx, col_idx-1, col_idx+1
        called_del <- FALSE

        # helper to test a candidate column
        check_gap_at <- function(j) {
          if (j < 1 || j > length(target_chars)) return(FALSE)
          # We only accept a gap if the reference residue at `pos` is indeed `ref_AA`
          # (i.e., the mutation spec matches the reference sequence)
          if (ref_chars[col_idx] != ref_AA) return(FALSE)
          # Count deletion if the read shows a gap at j and the reference at j is *not* a gap
          # (so it's a real residue position, not a ref gap)
          target_chars[j] == "-" && ref_chars[j] != "-"
        }

        if (check_gap_at(col_idx) ||
            check_gap_at(col_idx - 1) ||
            check_gap_at(col_idx + 1)) {
          score <- score + 1
          detected_mutations <- c(detected_mutations, paste0(ref_AA, pos, "-"))
          called_del <- TRUE
        }

        # If no deletion called but we do see explicit WT at the mapped column, confirm WT
        if (!called_del && aa_read_at_pos == ref_AA) {
          wt_confirmed <- wt_confirmed + 1
        }
        next
      }

      # --- Substitutions & WT confirmation ---
      # Treat gap/unknown in the READ as uninformative
      if (aa_read_at_pos %in% c("-", "X")) next

      if (aa_read_at_pos == target_AA) {
        # Resistance substitution detected
        score <- score + 1
        detected_mutations <- c(detected_mutations, paste0(ref_AA, pos, target_AA))
      } else if (aa_ref_at_pos == ref_AA && aa_read_at_pos == ref_AA) {
        # Explicit WT confirmation at this mutation site
        wt_confirmed <- wt_confirmed + 1
      }
      # else: some other mismatch – ignore under simplified rule
    }

    detected_mutations_str <- if (length(detected_mutations) > 0) paste(detected_mutations, collapse = ", ") else "None"

    # Status rule: Resistant if any R; else WT if any WT; else Unknown
    status <- if (score > 0) "Resistant" else if (wt_confirmed > 0) "Wildtype" else "Unknown"

    list(
      score = score,
      detected_mutations = detected_mutations_str,
      wt_confirmed_positions = wt_confirmed,
      status = status
    )

  }, error = function(e) {
    warning(paste("Error processing sequences:", e$message))
    list(score = NA, detected_mutations = "Error", wt_confirmed_positions = NA, status = "Error")
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

# Example usage:
result <- process_mutation_data("input.tsv")
print(result)
