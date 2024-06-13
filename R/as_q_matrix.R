#' Convert a Q matrix into a `q_matrix` obejct
#'
#' Takes a matrix of Q values, check its validity, and then formats it correctly
#' to make sure it can then be processed and plotted correctly
#'
#' @param x a matrix
#' @returns a `q_matrix` object, which is a matrix with appropriate column names
#' (.QX, where X is the component number) to use with plotting
#' @export

as_q_matrix <- function(x){
  if (inherits(x,"data.frame")){
    x <- as.matrix(x)
  }
  colnames(x)<- paste0(".Q",seq_len(ncol(x)))
  class(x) <- c("q_matrix",class(x))
  x
}





#' Tidy a Q matrix
#'
#' Takes a `q_matrix` object, which is a matrix, and returns a tidied tibble.
#'
#' @param x A Q matrix object (as returned by LEA::Q()).
#' @param gen_tbl An associated gen_tibble
#' @param ... not currently used
#' @return A tidied tibble
#' @export
tidy.q_matrix <- function(x, gen_tbl, ...){
  rlang::check_dots_empty()
  q_tbl <- x %>%
    tibble::as_tibble() %>%
    dplyr::mutate(id = gen_tbl$id,
                  group = gen_tbl$population)

  q_tbl <- q_tbl %>% tidyr::pivot_longer(cols = dplyr::starts_with(".Q"),
                                         names_to = "q", values_to = "percentage")
  #q_tbl
  dominant_q <- q_tbl %>%
    dplyr::group_by(.data$id) %>%
    dplyr::summarize(dominant_pop = .data$group[which.max(.data$percentage)], dominant_q = max(.data$percentage), .groups = 'drop')

  q_tbl <- q_tbl %>%
    dplyr::left_join(dominant_q, by = "id")

  q_tbl <- q_tbl %>%
    dplyr::arrange(.data$group, dplyr::desc(.data$dominant_q)) %>%
    dplyr::mutate(plot_order = dplyr::row_number(),  # Create plot_order column
           id = factor(.data$id, levels = unique(.data$id[order(.data$group, -.data$dominant_q)])))
}


#' Tidy ADMXITURE output files into plots
#'
#' Takes the name of a directory containing .Q file outputs, and
#' produces a list of tidied tibbles ready to plot.
#'
#' @param x the name of a directory containing .Q files
#' @param gen_tbl  An associated gen_tibble
#' @returns a list of `q_matrix` objects to plot
#'
#' @export

read_q_matrix_list <- function(x, gen_tbl){

  files <- list.files(x, pattern = "\\.Q$", full.names = TRUE)

  # Read all .Q files into a list
  data_list <- lapply(files, function(file) utils::read.table(file, header = FALSE))

  # Turn each into a Q matrix
  matrix_list <- lapply(data_list, FUN = as_q_matrix)

  # Sort matrix_list by the number of columns
  matrix_list <- matrix_list[order(sapply(matrix_list, ncol))]

  # Tidy each
  matrix_list <- lapply(matrix_list, function(x) tidy.q_matrix(x, gen_tbl = gen_tbl))

  matrix_list
}










