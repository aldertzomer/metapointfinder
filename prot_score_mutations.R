# prot_score_mutations.R — uses pwalign (global-local), unique WT counting, coverage guard

suppressPackageStartupMessages({
  library(Biostrings)
  library(pwalign)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
threads <- as.integer(args[2])


# -------- Core scorer (protein) --------
calculate_mutation_score <- function(reference, target, changes_str) {
  if (is.null(changes_str) || changes_str == "" || is.na(changes_str)) {
    return(list(score = 0, detected_mutations = "None", wt_confirmed_positions = 0, status = "Unknown"))
  }

  # Parse tokens like: A51V, Q42-, ASLD153-, R335RVPYR, N474KG, G194GSG
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

  # Normalize sequences
  reference <- toupper(gsub("[ \r\n\t]", "", reference))
  target    <- toupper(gsub("[- \r\n\t]", "", target))

  # Pairwise alignment (protein), full reference vs local subject
  data(BLOSUM62)
  pa <- pairwiseAlignment(
    pattern = AAString(reference),
    subject = AAString(target),
    type = "global-local",
    substitutionMatrix = BLOSUM62,
    gapOpening = -10, gapExtension = -0.5
  )

  aligned_ref    <- as.character(alignedPattern(pa))
  aligned_target <- as.character(alignedSubject(pa))
  ref_chars  <- strsplit(aligned_ref, "")[[1]]
  targ_chars <- strsplit(aligned_target, "")[[1]]

  # Map ungapped reference positions -> alignment columns
  ref_pos_to_col <- which(ref_chars != "-")
  max_ref_pos <- length(ref_pos_to_col)

  # Subject coverage window (columns where subject is not a gap)
  covered_cols <- which(targ_chars != "-")
  if (length(covered_cols) == 0L) {
    return(list(score = 0, detected_mutations = "None", wt_confirmed_positions = 0, status = "Unknown"))
  }
  cov_start <- min(covered_cols)
  cov_end   <- max(covered_cols)
  inside_cov <- function(cols) {
    length(cols) > 0L && min(cols) >= cov_start && max(cols) <= cov_end
  }

  take_at_cols   <- function(chars, cols) paste(chars[cols], collapse = "")
  is_informative <- function(aa) !(aa %in% c("-", "X"))

  score <- 0L
  detected <- character(0)

  # NEW: track unique positions for detected mutations and WT confirmations
  detected_pos <- integer(0)
  wt_pos <- integer(0)
  add_wt_once <- function(pos) {
    if (length(pos) == 1L && !(pos %in% detected_pos) && !(pos %in% wt_pos)) {
      wt_pos <<- c(wt_pos, pos)
    }
  }

  for (p in parsed) {
    pos     <- p$Pos
    ref_seg <- p$RefSeg
    alt_seg <- p$AltSeg
    raw_tok <- p$Raw

    ref_len <- nchar(ref_seg)
    alt_len <- nchar(alt_seg)

    # Position sanity
    if (is.na(pos) || pos < 1 || (pos + ref_len - 1) > max_ref_pos) next

    # Alignment columns spanned by the reference segment (contiguous in reference coords)
    ref_cols     <- ref_pos_to_col[pos:(pos + ref_len - 1)]
    ref_at_cols  <- take_at_cols(ref_chars, ref_cols)
    targ_at_cols <- take_at_cols(targ_chars, ref_cols)

    # Coverage guard: skip tokens outside the subject's aligned coverage
    if (!inside_cov(ref_cols)) next

    # ----- CASE 1: Deletions -----
    if (alt_seg == "-") {
      if (ref_len == 1) {
        col_idx <- ref_cols[1]
        if (ref_chars[col_idx] == ref_seg) {
          check_gap_at <- function(j) {
            if (j < 1 || j > length(targ_chars)) return(FALSE)
            (ref_chars[j] != "-") && (targ_chars[j] == "-")
          }
          if (check_gap_at(col_idx) || check_gap_at(col_idx - 1) || check_gap_at(col_idx + 1)) {
            score <- score + 1L
            detected <- c(detected, paste0(ref_seg, pos, "-"))
            detected_pos <- c(detected_pos, pos)
          } else if (targ_chars[col_idx] == ref_seg) {
            add_wt_once(pos)
          }
        }
      } else {
        if (ref_at_cols == ref_seg) {
          gaps <- vapply(ref_cols, function(j) targ_chars[j] == "-", logical(1))
          if (all(gaps)) {
            score <- score + 1L
            detected <- c(detected, paste0(ref_seg, pos, "-"))
            detected_pos <- c(detected_pos, pos)
          } else if (targ_at_cols == ref_seg) {
            add_wt_once(pos)
          }
        }
      }
      next
    }

    # ----- CASE 2: "Insertion-after" (e.g., R335RVPYR) -----
    if (ref_len == 1 && alt_len > 1 && substr(alt_seg, 1, 1) == ref_seg) {
      col_idx <- ref_cols[1]
      if (ref_chars[col_idx] == ref_seg) {
        base_ok <- (targ_chars[col_idx] == ref_seg)
        alt_tail <- substring(alt_seg, 2)
        tail_needed <- nchar(alt_tail)
        j <- col_idx + 1
        tail_seen <- character(0)
        while (j <= length(ref_chars) && nchar(paste(tail_seen, collapse = "")) < tail_needed) {
          if (ref_chars[j] == "-") {
            aa <- targ_chars[j]
            if (is_informative(aa)) tail_seen <- c(tail_seen, aa)
          } else {
            if (length(tail_seen) > 0) break
          }
          j <- j + 1
        }
        tail_seq <- paste(tail_seen, collapse = "")
        if (base_ok && tail_seq == alt_tail) {
          score <- score + 1L
          detected <- c(detected, raw_tok)
          detected_pos <- c(detected_pos, pos)
        } else if (base_ok) {
          add_wt_once(pos)
        }
      }
      next
    }

    # ----- CASE 3: Equal-length replacement/substitution (incl. multi-AA) -----
    if (ref_len == alt_len && alt_seg != ref_seg) {
      if (ref_at_cols == ref_seg) {
        targ_block <- targ_at_cols
        if (!grepl("[-X]", targ_block)) {
          if (targ_block == alt_seg) {
            score <- score + 1L
            detected <- c(detected, raw_tok)
            detected_pos <- c(detected_pos, pos)
          } else if (targ_block == ref_seg) {
            add_wt_once(pos)
          }
        }
      } else if (targ_at_cols == ref_seg && !grepl("[-X]", targ_at_cols)) {
        add_wt_once(pos)
      }
      next
    }

    # ----- CASE 4: Single-AA substitution fallback -----
    if (ref_len == 1 && alt_len == 1 && alt_seg != ref_seg) {
      col_idx <- ref_cols[1]
      aa_ref_at_pos  <- ref_chars[col_idx]
      aa_read_at_pos <- targ_chars[col_idx]
      if (is_informative(aa_read_at_pos)) {
        if (aa_ref_at_pos == ref_seg && aa_read_at_pos == alt_seg) {
          score <- score + 1L
          detected <- c(detected, paste0(ref_seg, pos, alt_seg))
          detected_pos <- c(detected_pos, pos)
        } else if (aa_ref_at_pos == ref_seg && aa_read_at_pos == ref_seg) {
          add_wt_once(pos)
        }
      }
      next
    }

    # ----- CASE 5: Other patterns not explicitly modeled -----
    if (ref_at_cols == ref_seg && targ_at_cols == ref_seg && !grepl("[-X]", targ_at_cols)) {
      add_wt_once(pos)
    }
  }

  detected_str <- if (length(detected) > 0) paste(detected, collapse = ", ") else "None"
  wt_conf <- length(unique(wt_pos))
  status <- if (score > 0) "Resistant" else if (wt_conf > 0) "Wildtype" else "Unknown"

  list(
    score = score,
    detected_mutations = detected_str,
    wt_confirmed_positions = wt_conf,
    status = status
  )
}

# -------- Runner: reads TSV and writes annotated output --------
process_mutation_data <- function(file_path) {
  input_data <- read.delim(file_path,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           quote = "",
                           check.names = FALSE)

  required_cols <- c("class", "gene", "read", "reference", "target", "changes_str")
  missing_cols <- setdiff(required_cols, colnames(input_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

#  results <- lapply(seq_len(nrow(input_data)), function(i) {
  results <- mclapply(seq_len(nrow(input_data)), function(i) {
    calculate_mutation_score(
      reference   = input_data$reference[i],
      target      = input_data$target[i],
      changes_str = input_data$changes_str[i]
    )
  }, mc.cores = threads)

  input_data$MutationScore         <- sapply(results, `[[`, "score")
  input_data$DetectedMutations     <- sapply(results, `[[`, "detected_mutations")
  input_data$WTConfirmedPositions  <- sapply(results, `[[`, "wt_confirmed_positions")
  input_data$Status                <- sapply(results, `[[`, "status")

  write.table(input_data,
              "updated_table_with_scores_and_mutations.tsv",
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)

  invisible(input_data)
}

# Example usage:
result <- process_mutation_data(input_file)
print(result)
