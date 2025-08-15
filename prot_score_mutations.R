# Required libraries
library(msa)
library(Biostrings)

# ----- Core scorer (protein) -----
# Reference-anchored indexing + robust indel/substitution handling with strict ref checks
calculate_mutation_score <- function(reference, target, changes_str) {
  if (is.null(changes_str) || changes_str == "" || is.na(changes_str)) {
    return(list(score = 0, detected_mutations = "None", wt_confirmed_positions = 0, status = "Unknown"))
  }

  # --- Parse tokens like A51V, Q42-, ASLD153-, R335RVPYR, N474KG, G194GSG ---
  tokens <- trimws(strsplit(changes_str, ",")[[1]])
  parse_token <- function(tok) {
    tok <- toupper(tok)
    m <- regexec("^([A-Z]+)([0-9]+)([A-Z\\-]+)$", tok)
    mm <- regmatches(tok, m)[[1]]
    if (length(mm) != 4) return(NULL)
    list(RefSeg = mm[2], Pos = as.integer(mm[3]), AltSeg = mm[4], Raw = tok)
  }
  parsed <- lapply(tokens, parse_token)
  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  if (length(parsed) == 0) {
    return(list(score = 0, detected_mutations = "None", wt_confirmed_positions = 0, status = "Unknown"))
  }

  tryCatch({
    # --- Normalize sequences (strip gaps from target); pad TARGET with X to avoid terminal artefacts ---
    reference <- toupper(gsub("[ \r\n\t]", "", reference))
    target    <- toupper(gsub("[- \r\n\t]", "", target))

    pad_len <- 2000L
    pad_blk <- paste(rep("X", pad_len), collapse = "")
    target  <- paste0(pad_blk, target, pad_blk)

    # --- Align (ClustalW keeps flanks) ---
    sequences <- AAStringSet(c(Reference = reference, Target = target))
    aln <- msa(sequences, method = "ClustalW", verbose = FALSE)
    aligned_ref    <- as.character(unmasked(aln)["Reference"])
    aligned_target <- as.character(unmasked(aln)["Target"])
    ref_chars    <- strsplit(aligned_ref, "")[[1]]
    targ_chars   <- strsplit(aligned_target, "")[[1]]

    # Map ungapped reference positions -> alignment columns
    ref_pos_to_col <- which(ref_chars != "-")
    max_ref_pos <- length(ref_pos_to_col)

    # Helpers
    take_at_cols   <- function(chars, cols) paste(chars[cols], collapse = "")
    is_informative <- function(aa) !(aa %in% c("-", "X"))

    score <- 0L
    detected <- character(0)
    wt_conf <- 0L

    for (p in parsed) {
      pos     <- p$Pos
      ref_seg <- p$RefSeg
      alt_seg <- p$AltSeg
      raw_tok <- p$Raw

      ref_len <- nchar(ref_seg)
      alt_len <- nchar(alt_seg)

      if (is.na(pos) || pos < 1 || (pos + ref_len - 1) > max_ref_pos) next

      # Alignment columns spanned by the reference segment (contiguous in reference coords)
      ref_cols     <- ref_pos_to_col[pos:(pos + ref_len - 1)]
      ref_at_cols  <- take_at_cols(ref_chars, ref_cols)
      targ_at_cols <- take_at_cols(targ_chars, ref_cols)

      # ========== CASE 1: Deletions ==========
      if (alt_seg == "-") {
        if (ref_len == 1) {
          # Single-AA deletion, allow ±1 col jitter but REQUIRE correct ref at site
          col_idx <- ref_cols[1]
          # require that the reference column holds the expected ref AA
          if (ref_chars[col_idx] == ref_seg) {
            check_gap_at <- function(j) {
              if (j < 1 || j > length(targ_chars)) return(FALSE)
              (ref_chars[j] != "-") && (targ_chars[j] == "-")
            }
            if (check_gap_at(col_idx) || check_gap_at(col_idx - 1) || check_gap_at(col_idx + 1)) {
              score <- score + 1L
              detected <- c(detected, paste0(ref_seg, pos, "-"))
            } else if (targ_chars[col_idx] == ref_seg) {
              wt_conf <- wt_conf + 1L
            }
          }
        } else {
          # Multi-AA deletion: require token's ref segment is present at those ref columns
          # and all those columns are gaps in the target
          if (ref_at_cols == ref_seg) {
            gaps <- vapply(ref_cols, function(j) targ_chars[j] == "-", logical(1))
            if (all(gaps)) {
              score <- score + 1L
              detected <- c(detected, paste0(ref_seg, pos, "-"))
            } else if (targ_at_cols == ref_seg) {
              wt_conf <- wt_conf + 1L
            }
          }
        }
        next
      }

      # ========== CASE 2: "Insertion-after" (e.g., R335RVPYR) ==========
      # Keep first ref AA (must match), insert tail immediately after in reference-gap columns.
      if (ref_len == 1 && alt_len > 1 && substr(alt_seg, 1, 1) == ref_seg) {
        col_idx <- ref_cols[1]
        # Strict: reference at this column must match the token's ref AA
        if (ref_chars[col_idx] == ref_seg) {
          base_ok <- (targ_chars[col_idx] == ref_seg)
          alt_tail <- substring(alt_seg, 2)
          tail_needed <- nchar(alt_tail)
          # scan forward to collect insertion residues only where reference has '-'
          j <- col_idx + 1
          tail_seen <- character(0)
          while (j <= length(ref_chars) && nchar(paste(tail_seen, collapse = "")) < tail_needed) {
            if (ref_chars[j] == "-") {
              aa <- targ_chars[j]
              if (is_informative(aa)) tail_seen <- c(tail_seen, aa)
            } else {
              # once we re-enter real reference, stop if we've started collecting
              if (length(tail_seen) > 0) break
            }
            j <- j + 1
          }
          tail_seq <- paste(tail_seen, collapse = "")
          if (base_ok && tail_seq == alt_tail) {
            score <- score + 1L
            detected <- c(detected, raw_tok)
          } else if (base_ok) {
            wt_conf <- wt_conf + 1L
          }
        }
        next
      }

      # ========== CASE 3: Equal-length replacement/substitution (incl. multi-AA) ==========
      if (ref_len == alt_len && alt_seg != ref_seg) {
        # STRICT: require the reference block equals the token's RefSeg
        if (ref_at_cols == ref_seg) {
          targ_block <- targ_at_cols
          # Need all informative residues to compare
          if (!grepl("[-X]", targ_block)) {
            if (targ_block == alt_seg) {
              score <- score + 1L
              detected <- c(detected, raw_tok)
            } else if (targ_block == ref_seg) {
              wt_conf <- wt_conf + 1L
            }
          }
        } else if (targ_at_cols == ref_seg && !grepl("[-X]", targ_at_cols)) {
          # explicit WT block at site
          wt_conf <- wt_conf + 1L
        }
        next
      }

      # ========== CASE 4: Single-AA substitution fallback with strict ref check ==========
      if (ref_len == 1 && alt_len == 1 && alt_seg != ref_seg) {
        col_idx <- ref_cols[1]
        aa_ref_at_pos  <- ref_chars[col_idx]
        aa_read_at_pos <- targ_chars[col_idx]
        # Treat gap/unknown in READ as uninformative
        if (!(aa_read_at_pos %in% c("-", "X"))) {
          # STRICT: ref at site must match token's ref AA
          if (aa_ref_at_pos == ref_seg && aa_read_at_pos == alt_seg) {
            score <- score + 1L
            detected <- c(detected, paste0(ref_seg, pos, alt_seg))
          } else if (aa_ref_at_pos == ref_seg && aa_read_at_pos == ref_seg) {
            wt_conf <- wt_conf + 1L
          }
        }
        next
      }

      # ========== CASE 5: Other patterns not explicitly modeled ==========
      # Be conservative: only confirm WT if the exact ref segment is present and informative
      if (ref_at_cols == ref_seg && targ_at_cols == ref_seg && !grepl("[-X]", targ_at_cols)) {
        wt_conf <- wt_conf + 1L
      }
      # else: leave as Unknown
    }

    detected_str <- if (length(detected) > 0) paste(detected, collapse = ", ") else "None"
    status <- if (score > 0) "Resistant" else if (wt_conf > 0) "Wildtype" else "Unknown"

    list(
      score = score,
      detected_mutations = detected_str,
      wt_confirmed_positions = wt_conf,
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
