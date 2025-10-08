# Required libraries
library(Biostrings)
library(pwalign)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]         # e.g. "<sample>.prot.input.tsv"
threads    <- as.integer(args[2])

calculate_mutation_score <- function(reference, target, changes_str, debug = FALSE) {
  if (is.null(changes_str) || changes_str == "" || is.na(changes_str)) {
    return(list(score = 0, detected_mutations = "None",
                wt_confirmed_positions = 0, status = "Unknown"))
  }

  # Parse tokens like: C12A, T22G, G39GG, T40TGT, C43T, G43-, etc.
  change_tokens <- trimws(strsplit(changes_str, ",")[[1]])
  parse_token <- function(tok) {
    tok <- toupper(tok)
    m <- regexec("^([ACGT]+)([0-9]+)([ACGT\\-]+)$", tok)
    mm <- regmatches(tok, m)[[1]]
    if (length(mm) != 4) return(NULL)
    list(RefSeg = mm[2], Pos = as.integer(mm[3]), AltSeg = mm[4], Raw = tok)
  }
  parsed <- lapply(change_tokens, parse_token)
  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  if (length(parsed) == 0) {
    return(list(score = 0, detected_mutations = "None",
                wt_confirmed_positions = 0, status = "Unknown"))
  }

  # Normalize inputs
  reference <- toupper(gsub("[ \r\n\t]", "", reference))
  target    <- toupper(gsub("[- \r\n\t]", "", target))

  # Align: force full reference (pattern) to align; read (subject) can be local
  submat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -5, baseOnly = FALSE)  # tolerate N/IUPAC
  pa <- pairwiseAlignment(
    pattern = DNAString(reference),
    subject = DNAString(target),
    type = "global-local",
    substitutionMatrix = submat,
    gapOpening = -10, gapExtension = -0.5
  )

  aligned_ref    <- as.character(alignedPattern(pa))
  aligned_target <- as.character(alignedSubject(pa))
  ref_chars    <- strsplit(aligned_ref, "")[[1]]
  target_chars <- strsplit(aligned_target, "")[[1]]

  # Build map: alignment column -> (ungapped) reference position
  ref_pos_at_col <- rep(NA_integer_, length(ref_chars))
  rp <- 0L
  for (k in seq_along(ref_chars)) {
    if (ref_chars[k] != "-") {
      rp <- rp + 1L
      ref_pos_at_col[k] <- rp
    }
  }
  max_ref_pos <- max(ref_pos_at_col, na.rm = TRUE)

  informative <- function(b) !(b %in% c("-", "N"))

  # Subject coverage window (columns where subject is not a gap)
  covered_cols <- which(target_chars != "-")
  if (length(covered_cols) == 0L) {
    return(list(score = 0, detected_mutations = "None",
                wt_confirmed_positions = 0, status = "Unknown"))
  }
  cov_start <- min(covered_cols); cov_end <- max(covered_cols)
  inside_cov <- function(cols_or_j) {
    if (length(cols_or_j) == 1L) {
      j <- cols_or_j
      !is.na(j) && j >= cov_start && j <= cov_end
    } else {
      length(cols_or_j) > 0L && min(cols_or_j) >= cov_start && max(cols_or_j) <= cov_end
    }
  }

  take_at_cols <- function(chars, cols) paste(chars[cols], collapse = "")

  score <- 0L
  detected <- character(0)

  # Track unique positions for WT confirmation and detected calls
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

    # Sanity on reference coordinate range
    if (is.na(pos) || pos < 1 || (pos + ref_len - 1) > max_ref_pos) next

    # Alignment columns spanned by the reference segment
    ref_cols <- which(ref_pos_at_col %in% pos:(pos + ref_len - 1))
    if (length(ref_cols) != ref_len) next  # ambiguous/misaligned
    if (!inside_cov(ref_cols)) next        # coverage guard

    ref_block <- take_at_cols(ref_chars, ref_cols)
    read_block <- take_at_cols(target_chars, ref_cols)

    # Reference spec must match alignment's reference block
    if (ref_block != ref_seg) next

    # ---- CASE A: Deletion ----
    if (alt_seg == "-") {
      gaps <- vapply(ref_cols, function(j) target_chars[j] == "-", logical(1))
      if (all(gaps)) {
        score <- score + 1L
        detected <- c(detected, paste0(ref_seg, pos, "-"))
        detected_pos <- c(detected_pos, pos)
      } else if (read_block == ref_seg && !grepl("N", read_block)) {
        add_wt_once(pos)
      }
      next
    }

    # ---- CASE B: Insertion-after (e.g., G39GG, T40TGT) ----
    if (ref_len == 1 && alt_len > 1 && startsWith(alt_seg, ref_seg)) {
      j <- ref_cols[1]
      base_ok <- (target_chars[j] == ref_seg)
      alt_tail <- substring(alt_seg, 2)            # bases inserted after the ref base
      tail_needed <- nchar(alt_tail)

      # Walk right collecting subject bases where reference has gaps
      jj <- j + 1
      tail_seen <- character(0)
      started <- FALSE
      while (jj <= length(ref_chars) && nchar(paste0(tail_seen, collapse = "")) < tail_needed) {
        if (ref_chars[jj] == "-") {
          aa <- target_chars[jj]
          if (informative(aa)) {
            tail_seen <- c(tail_seen, aa)
            started <- TRUE
          }
        } else {
          if (started) break  # stop at first ref letter after starting insertion tract
        }
        jj <- jj + 1
      }
      tail_seq <- paste0(tail_seen, collapse = "")
      if (base_ok && tail_seq == alt_tail) {
        score <- score + 1L
        detected <- c(detected, raw_tok)
        detected_pos <- c(detected_pos, pos)
      } else if (base_ok) {
        add_wt_once(pos)
      }
      next
    }

    # ---- CASE C: Equal-length replacement (incl. multi-base) ----
    if (ref_len == alt_len && alt_seg != ref_seg) {
      if (!grepl("[-N]", read_block) && read_block == alt_seg) {
        score <- score + 1L
        detected <- c(detected, raw_tok)
        detected_pos <- c(detected_pos, pos)
      } else if (!grepl("[-N]", read_block) && read_block == ref_seg) {
        add_wt_once(pos)
      }
      next
    }

    # ---- CASE D: Single-base substitution fallback ----
    if (ref_len == 1 && alt_len == 1 && alt_seg != ref_seg) {
      j <- ref_cols[1]
      rj <- ref_chars[j]
      sj <- target_chars[j]
      if (informative(sj)) {
        if (rj == ref_seg && sj == alt_seg) {
          score <- score + 1L
          detected <- c(detected, paste0(ref_seg, pos, alt_seg))
          detected_pos <- c(detected_pos, pos)
        } else if (rj == ref_seg && sj == ref_seg) {
          add_wt_once(pos)
        }
      }
      next
    }

    # ---- CASE E: Nothing matched → possible WT confirmation ----
    if (!grepl("[-N]", read_block) && read_block == ref_seg) {
      add_wt_once(pos)
    }
  }

  detected_str <- if (length(detected) > 0) paste(detected, collapse = ",") else "None"
  wt_confirmed <- length(unique(wt_pos))
  status <- if (score > 0) "Resistant" else if (wt_confirmed > 0) "Wildtype" else "Unknown"

  list(score = score,
       detected_mutations = detected_str,
       wt_confirmed_positions = wt_confirmed,
       status = status)
}

process_mutation_data <- function(file_path) {
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

  input_data
}

# Run with the same workflow as before
result <- process_mutation_data(input_file)
warnings()
print(result)
