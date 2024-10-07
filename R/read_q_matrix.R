read_q_matrix <- function(x) {
  # Check if the input is already a read-in data frame or matrix
  if (inherits(x, "data.frame") || inherits(x, "matrix")) {
    message("Input is a data frame or matrix")
    x <- as.matrix(x)
    colnames(x) <- paste0(".Q", seq_len(ncol(x)))
    class(x) <- c("q_matrix", class(x))
    return(x)
  }
  # Check if the input exists as a file or directory
  if (file.exists(x)) {
    if (dir.exists(x)) {
      message("Input is a directory")
      # List all .Q files in the directory
      files <- list.files(x, pattern = "\\.Q$", full.names = TRUE)
      # Check if the directory contains at least one .Q file
      if (length(files) == 0) {
        stop("No .Q files found in the directory")
      }
      # Read all .Q files into a list
      data_list <- lapply(files, function(file) utils::read.table(file, header = FALSE))
      # Turn each data frame into a Q matrix
      matrix_list <- lapply(data_list, FUN = as_q_matrix)
      # Get the number of columns for each Q matrix (which corresponds to K)
      k_values <- sapply(matrix_list, ncol)
      # Create a list of lists, grouping matrices by their number of columns (K)
      list_of_lists <- split(matrix_list, k_values)
      # Rename the list so that each K is prefixed with "k"
      names(list_of_lists) <- paste0("k", names(list_of_lists))
      # Return the list of lists
     list_of_lists
    } else {
      # Action if input is a file
      message("Input is a file")
      if (grepl("\\.Q$", x)) {
        x <- utils::read.table(x, header = FALSE)
        x <- as.matrix(x)
        colnames(x) <- paste0(".Q", seq_len(ncol(x)))
        class(x) <- c("q_matrix", class(x))
        return(x)
      } else {
        #if input file does not end in Q
        stop("Input file does not end in '.Q'")
      }
    }
  } else {
    # If the input is not a valid file or directory
    stop("Input is not a valid dataframe, file or directory")
  }
}



# Function to summarize the number of different K values and repeats
summarise_q_matrices <- function(list_of_lists) {
  # Get the K values (names of the list_of_lists)
  k_values <- names(list_of_lists)

  # Initialize an empty data frame to store the summary
  summary_df <- data.frame(K = integer(), Repeats = integer())

  # Loop through each K value and count the number of repeats (length of the sublist)
  for (k in k_values) {
    # Remove the "k" prefix and convert to numeric
    k_numeric <- as.numeric(sub("k", "", k))

    # Count the number of repeats for this K
    num_repeats <- length(list_of_lists[[k]])

    # Add the K and the number of repeats to the summary data frame
    summary_df <- rbind(summary_df, data.frame(K = k_numeric, Repeats = num_repeats))
  }

  return(summary_df)
}

















