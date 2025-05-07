# Required libraries
library(msa)
library(Biostrings)

# Function to calculate mutation score and detect mutations
# (Previous function remains the same)
calculate_mutation_score <- function(reference, target, changes_str) {
  # Handle empty or invalid input
  if (is.null(changes_str) || changes_str == "" || is.na(changes_str)) {
    return(list(score = 0, detected_mutations = "None"))
  }
  
  # Parse the comma-delimited changes string into a list
  change_pairs <- strsplit(changes_str, ",")[[1]]
  change_pairs <- trimws(change_pairs)  # Remove any whitespace
  
  # Create a data frame to store position, reference base, and target base
  change_data <- data.frame(
    Position = numeric(0),
    ReferenceDNA = character(0),
    TargetDNA = character(0),
    stringsAsFactors = FALSE
  )
  
  # Populate the change_data data frame
  for (change in change_pairs) {
    # Extract position and reference/target bases from each change (e.g., A51V -> Position 51, A -> reference, V -> target)
    # Add error handling for malformed change strings
    if (nchar(change) < 3) next
    
    position <- as.numeric(gsub("[^0-9]", "", change))
    if (is.na(position)) next
    
    reference_DNA <- substr(change, 1, 1)
    target_DNA <- substr(change, nchar(change), nchar(change))
    
    # Add to the data frame
    change_data <- rbind(change_data, 
                        data.frame(Position = position, 
                                 ReferenceDNA = reference_DNA, 
                                 TargetDNA = target_DNA))
  }
  
  # Convert sequences to DNAStringSet objects with proper error handling
  tryCatch({
    sequences <- DNAStringSet(c(Reference = reference, Target = target))
    alignment <- msa(sequences, method = "ClustalW", verbose = FALSE)
    
    # Extract aligned sequences
    aligned_ref <- as.character(unmasked(alignment)["Reference"])
    aligned_target <- as.character(unmasked(alignment)["Target"])
    
    # Split into character vectors
    ref_chars <- strsplit(aligned_ref, "")[[1]]
    target_chars <- strsplit(aligned_target, "")[[1]]
    
    # Initialize score and detected mutations vector
    score <- 0
    detected_mutations <- character(0)
    
    # Iterate through the change data and check for mutations in the target sequence
    for (i in 1:nrow(change_data)) {
      pos <- change_data$Position[i]
      ref_DNA <- change_data$ReferenceDNA[i]
      target_DNA <- change_data$TargetDNA[i]
      
      # Ensure position is within bounds of the target sequence
      if (pos <= length(target_chars)) {
        # Check if the target base at the position matches the expected change
        if (target_chars[pos] == target_DNA) {
          score <- score + 1
          detected_mutations <- c(detected_mutations, 
                                paste0(ref_DNA, pos, target_DNA))
        }
      }
    }
    
    # Format detected mutations string
    detected_mutations_str <- if (length(detected_mutations) > 0) {
      paste(detected_mutations, collapse = ", ")
    } else {
      "None"
    }
    
    return(list(
      score = score,
      detected_mutations = detected_mutations_str
    ))
  }, error = function(e) {
    warning(paste("Error processing sequences:", e$message))
    return(list(
      score = NA,
      detected_mutations = "Error"
    ))
  })
}

# Updated function to read and process the input file
process_mutation_data <- function(file_path) {
  # Read the tab-delimited file with error handling
  tryCatch({
    input_data <- read.delim(file_path, 
                            sep = "\t",           # Tab delimiter
                            stringsAsFactors = FALSE,
                            quote = "",           # No quotes
                            check.names = FALSE)  # Keep column names as is
    
    # Verify required columns exist
    required_cols <- c("class", "gene", "read", "reference", "target", "changes_str")
    missing_cols <- required_cols[!required_cols %in% colnames(input_data)]
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    }
    
    # Calculate mutation scores and detect mutations
    results <- lapply(1:nrow(input_data), function(i) {
      calculate_mutation_score(
        input_data$reference[i],
        input_data$target[i],
        input_data$changes_str[i]
      )
    })
    
    # Extract scores and detected mutations
    mutation_scores <- sapply(results, function(x) x$score)
    detected_mutations <- sapply(results, function(x) x$detected_mutations)
    
    # Add results to the data frame
    input_data$MutationScore <- mutation_scores
    input_data$DetectedMutations <- detected_mutations
    
    # Write results (tab-delimited)
    output_file <- "updated_table_with_scores_and_mutations.tsv"
    write.table(input_data, 
                output_file, 
                sep = "\t",           # Tab delimiter
                row.names = FALSE,    # No row names
                quote = FALSE)        # No quotes
    
    return(input_data)
  }, error = function(e) {
    stop(paste("Error processing file:", e$message))
  })
}

# Example usage:
result <- process_mutation_data("input.tsv")
print(result)