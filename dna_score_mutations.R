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
    Position     = numeric(0),
    ReferenceDNA = character(0),
    TargetDNA    = character(0),
    stringsAsFactors = FALSE
  )

  for (chg in change_pairs) {
    chg <- toupper(chg)
    if (!grepl("^[ACGT]+[0-9]+[ACGT-]+$", chg)) next
    pos <- as.integer(gsub("^[ACGT]+([0-9]+)[ACGT-]+$", "\\1", chg))
    ref_base <- gsub("^([ACGT]+)[0-9]+[ACGT-]+$", "\\1", chg)
    if (nchar(ref_base) != 1) next
    ref_base <- substr(ref_base, 1, 1)
    tgt_base <- toupper(substr(chg, nchar(chg), nchar(chg)))  # "-" for deletion
    change_data <- rbind(change_data,
                         data.frame(Position = pos, ReferenceDNA = ref_base, TargetDNA = tgt_base))
  }

  tryCatch({
    # Normalize inputs
    reference <- toupper(gsub("[ \r\n\t]", "", reference))
    target    <- toupper(gsub("[- \r\n\t]", "", target))

    # Align: force full reference (pattern) to align; read (subject) can overhang
    submat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -5, baseOnly = FALSE)
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

    # Build alignment-column -> reference-position map (NA where ref has a gap)
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

    # --- NEW: subject coverage window
    covered_cols <- which(target_chars != "-")
    if (length(covered_cols) == 0L) {
      return(list(score = 0, detected_mutations = "None",
                  wt_confirmed_positions = 0, status = "Unknown"))
    }
    cov_start <- min(covered_cols); cov_end <- max(covered_cols)
    inside_cov <- function(j) !is.na(j) && j >= cov_start && j <= cov_end

    score <- 0L
    detected <- character(0)

    # --- NEW: track unique positions (WT counted once per site)
    detected_pos <- integer(0)
    wt_pos <- integer(0)
    add_wt_once <- function(pos) {
      if (length(pos) == 1L && !(pos %in% detected_pos) && !(pos %in% wt_pos)) {
        wt_pos <<- c(wt_pos, pos)
      }
    }

    # Evaluate each change
    for (i in seq_len(nrow(change_data))) {
      pos      <- change_data$Position[i]
      ref_DNA  <- change_data$ReferenceDNA[i]
      target_DNA <- change_data$TargetDNA[i]

      if (is.na(pos) || pos < 1 || pos > max_ref_pos) next

      # Find alignment column j that maps to this reference position
      j <- which(ref_pos_at_col == pos)
      if (length(j) != 1L) next  # ambiguous or missing
      j <- j[[1]]

      # --- NEW: coverage guard
      if (!inside_cov(j)) next

      ref_at_j  <- ref_chars[j]
      read_at_j <- target_chars[j]

      # sanity: spec must match reference
      if (ref_at_j != ref_DNA) {
        if (debug) message(sprintf("Spec/ref mismatch at pos %d: spec=%s ref=%s", pos, ref_DNA, ref_at_j))
        next
      }

      if (target_DNA == "-") {
        # Deletion: require a gap in subject exactly at the ref-mapped column
        if (read_at_j == "-" && ref_at_j != "-") {
          score <- score + 1L
          detected <- c(detected, paste0(ref_DNA, pos, "-"))
          detected_pos <- c(detected_pos, pos)
          if (debug) message(sprintf("DEL %s at pos %d (col %d)", paste0(ref_DNA, pos, "-"), pos, j))
        } else if (informative(read_at_j) && read_at_j == ref_DNA) {
          add_wt_once(pos)
          if (debug) message(sprintf("WT confirm (no deletion) at pos %d", pos))
        } else if (debug) {
          message(sprintf("No DEL at pos %d: read(col)=%s", pos, read_at_j))
        }
      } else {
        # Substitution: subject base must equal target_DNA at mapped column
        if (informative(read_at_j) && read_at_j == target_DNA) {
          score <- score + 1L
          detected <- c(detected, paste0(ref_DNA, pos, target_DNA))
          detected_pos <- c(detected_pos, pos)
          if (debug) message(sprintf("SUB %s at pos %d (col %d)", paste0(ref_DNA, pos, target_DNA), pos, j))
        } else if (informative(read_at_j) && read_at_j == ref_DNA) {
          add_wt_once(pos)
          if (debug) message(sprintf("WT confirm at pos %d", pos))
        } else if (debug) {
          message(sprintf("No SUB at pos %d: read(col)=%s", pos, read_at_j))
        }
      }
    }

    detected_str <- if (length(detected) > 0) paste(detected, collapse = ",") else "None"
    wt_confirmed <- length(unique(wt_pos))  # NEW: unique per position

    status <- if (score > 0) "Resistant" else if (wt_confirmed > 0) "Wildtype" else "Unknown"

    list(score = score,
         detected_mutations = detected_str,
         wt_confirmed_positions = wt_confirmed,
         status = status)

  }, error = function(e) {
    if (debug) message(paste("Error processing sequences:", e$message))
    return(list(score = NA, detected_mutations = "Error",
                wt_confirmed_positions = NA, status = "Error"))
  })
}

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

# Run
result <- process_mutation_data("input.tsv")
warnings()
print(result)
