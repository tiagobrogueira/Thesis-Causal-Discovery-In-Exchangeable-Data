library(tidyverse) # For data manipulation (read_csv, tibble, etc.)
library(Matrix)    # For sparse matrices if needed, though not directly used in the Python example
library(here)      # For robust file path management



run_synthetic <- function(func, dataset, func_name, ...) {
  # Define the directory containing the dataset
  dir <- paste0("generate_txt/", dataset)
  
  # Read all files in the directory
  files <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)
  
  # Initialize variables to store scores and weights
  scores <- numeric(length(files))
  weights <- rep(1 / length(files), length(files))  # Equal weights
  
  # Loop through each file
  for (i in seq_along(files)) {
    # Read the data (assuming 2 columns)
    data <- as.matrix(read.table(files[i]))
    if (any(is.na(data))) {
      print(paste("Found NaN in data:", dataset, files[i]))
      next
    }
    if (ncol(data) != 2) {
    message("Skipping file ", files[i], ": does not have exactly 2 columns.")
    scores[i] <- 0
    next
    }
  
    
    # Call the function and store the score
    res <- func(data, ...)
    scores[i] <- res$eps  # Assuming the score is in res$eps
  }
  
  # Create a folder for predictions
  kwargs_str <- paste(names(list(...)), unlist(list(...)), sep = "=", collapse = ", ")
  folder <- paste0("predictions/", func_name, "R_", kwargs_str)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  print("folder created")
  print(folder)
  print("func")
  print(func)
  
  # Save the scores and weights to a file
  output_file <- paste0(folder, "/", dataset, ".txt")
  write.table(cbind(scores, weights), file = output_file, row.names = FALSE, col.names = FALSE)
}

run_tuebingen <- function(func, ...) {
  # Define file paths and prefixes
  pairmeta_file <- "../Pairs/pairmeta.txt"
  pair_prefix <- "pair"
  
  # Read metadata
  meta_lines <- readLines(pairmeta_file)
  
  data_list <- list()
  weights <- numeric(length(meta_lines))
  
  for (i in seq_along(meta_lines)) {
    entries <- strsplit(meta_lines[i], "\\s+")[[1]]
    
    pair_number <- sprintf("%04d", as.integer(entries[1]))
    
    x_start <- as.integer(entries[2]) - 1
    x_end   <- as.integer(entries[3])
    y_start <- as.integer(entries[4]) - 1
    y_end   <- as.integer(entries[5])
    weight  <- as.numeric(entries[6])
    
    pair_file <- paste0("../Pairs/", pair_prefix, pair_number, ".txt")
    
    if (!file.exists(pair_file)) {
      message("Warning: ", pair_file, " not found, skipping.")
      next
    }
    
    pair_data <- tryCatch({
      as.matrix(read.table(pair_file))
    }, error = function(e) {
      message("Error reading ", pair_file, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(pair_data)) next
    
    # Extract submatrices
    x_data <- pair_data[, (x_start + 1):x_end, drop = FALSE]
    y_data <- pair_data[, (y_start + 1):y_end, drop = FALSE]
    
    data_list[[i]] <- list(x = x_data, y = y_data)
    weights[i] <- weight
  }
  
  # Remove NULL entries if any
  valid_indices <- which(!sapply(data_list, is.null))
  data_list <- data_list[valid_indices]
  weights <- weights[valid_indices]
  
  # Initialize scores
  scores <- numeric(length(data_list))
  
  # Compute scores using the provided function
  for (i in seq_along(data_list)) {
    xy_data <- cbind(data_list[[i]]$x, data_list[[i]]$y)
    
    if (any(is.na(xy_data))) {
      message("Found NaN in data for pair ", i)
      scores[i] <- 0
      next
    }
    
    res <- func(xy_data, ...)
    scores[i] <- res$eps
  }
  
  # Create folder and save results
  kwargs_str <- paste(names(list(...)), unlist(list(...)), sep = "=", collapse = ", ")
  folder <- paste0("predictions/", deparse(substitute(func)), "R_", kwargs_str)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  
  output_file <- paste0(folder, "/tuebingen.txt")
  write.table(cbind(scores, weights), file = output_file, row.names = FALSE, col.names = FALSE)
}

run_old <- function(dataset, func, ...) {
  # Define the folder path using the 'here' package for robustness
  # This assumes your R project root is where 'other implementations' is relative to.
  folder <- here("..", "other implementations", "synthetic_datasets")

  pairs_file <- file.path(folder, paste0(dataset, "_pairs.csv"))
  targets_file <- file.path(folder, paste0(dataset, "_targets.csv"))

  # --- sanity checks ---
  if (!file.exists(pairs_file)) {
    stop(paste("Pairs file not found:", pairs_file))
  }
  if (!file.exists(targets_file)) {
    stop(paste("Targets file not found:", targets_file))
  }

  # --- read in with pandas (using read_csv from readr package, part of tidyverse) ---
  df_pairs <- read_csv(pairs_file, show_col_types = FALSE)
  df_targets <- read_csv(targets_file, show_col_types = FALSE)

  if (nrow(df_pairs) != nrow(df_targets)) {
    stop(paste0("Row count mismatch: ", nrow(df_pairs), " in pairs vs ",
                nrow(df_targets), " in targets"))
  }

  data_list <- list()
  for (idx in seq_len(nrow(df_pairs))) {
    pair_row <- df_pairs[idx, ]
    target_row <- df_targets[idx, ]

    # get the raw space-separated strings
    # Assuming the second and third columns are the string representations,
    # similar to `iloc[1]` and `iloc[2]` in pandas (0-indexed in Python, 1-indexed in R)
    x_str <- as.character(pair_row[[2]]) # Use [[2]] for the second column
    y_str <- as.character(pair_row[[3]]) # Use [[3]] for the third column

    # split & convert to floats, then make column vectors (matrices in R)
    # Using as.numeric(strsplit(...)[[1]]) for robust conversion
    x <- matrix(as.numeric(unlist(strsplit(x_str, " "))), ncol = 1)
    y <- matrix(as.numeric(unlist(strsplit(y_str, " "))), ncol = 1)

    # swap if target's 2nd column is -1
    if (target_row[[2]] == -1) { # Using [[2]] for the second column of targets
      temp <- x
      x <- y
      y <- temp
    }
    # print(x) # Uncomment for debugging if needed
    # --- CORRECTION APPLIED HERE: Store x and y as named elements ---
    data_list[[idx]] <- list(x = x, y = y)
  }

  # `data` is already `data_list` here
  # Generate weights
  weights <- rep(1 / length(data_list), length(data_list))

  # Extract kwargs for folder naming
  kwargs <- list(...)
  kwargs_str <- paste(names(kwargs), unlist(kwargs), sep = "=", collapse = ", ")
  if (nchar(kwargs_str) > 0) {
    kwargs_str <- paste0("_", kwargs_str) # Prepend underscore if kwargs exist
  }

  # Construct folder name
  # func.__name__ in Python is equivalent to deparse(substitute(func)) in R
  output_folder <- paste0("predictions/", deparse(substitute(func)), "R_", kwargs_str)

  # Create directory if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

    scores <- numeric(length(data_list))
  
  # Compute scores using the provided function
  for (i in seq_along(data_list)) {
    xy_data <- cbind(data_list[[i]]$x, data_list[[i]]$y)
    
    if (any(is.na(xy_data))) {
      message("Found NaN in data for pair ", i)
      scores[i] <- 0
      next
    }
    
    res <- func(xy_data, ...)
    scores[i] <- res$eps
  }
  
  # Save scores and weights
  output_file_path <- file.path(output_folder, paste0(dataset, ".txt"))
  write.table(cbind(scores, weights),
              file = output_file_path, # Ensure this is output_file_path
              row.names = FALSE, col.names = FALSE)

}


test_all <- function(method) {
  # Path to the parent folder
  parent_folder <- "generate_txt/"
  
  # List all folders inside the parent folder
  datasets <- list.dirs(parent_folder, full.names = FALSE, recursive = FALSE)
  func_name <- deparse(substitute(method))
  starter = "exchangeable_log_raleigh_uni_"
  # Loop over each dataset folder
  for (dataset in datasets) {
    if (dataset < starter) {
      print("skipping")
      next
    }
    run_synthetic(method, dataset, func_name)
  }
}

source("../other implementations/bqcd/bqcd.R")
source("../other implementations/slope-20181208/Slope.R")
source("../other implementations/slope-20181208/utilities.R")

# Define a wrapper function for QCD_wrap
qcd_function <- function(data, ...) {
  QCD_wrap(data[,1], data[,2], "QCD_qnn")
}

Slope_ <- function(...) {
  res <- Slope(...)
  res$eps <- - res$eps
  return (res)
}

current_path <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(current_path))

 
run_tuebingen(qcd_function)
#test_all(qcd_function)
#run_old("CE-Net",qcd_function)
#run_old("CE-Net",Slope_)
#run_tuebingen(Slope_)
#test_all(Slope)
print("done")

