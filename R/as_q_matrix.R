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
#' @param data An associated tibble (e.g. a [`gen_tibble`]), with the individuals in the same order as the data used to
#' generate the Q matrix
#' @param ... not currently used
#' @return A tidied tibble
#' @export
tidy.q_matrix <- function(x, data, ...){
  rlang::check_dots_empty()
  q_tbl <- x %>%
    tibble::as_tibble() %>%
    dplyr::mutate(id = data$id,
                  # @TODO we should get this from the grouped tibble, not hardcode it!
                  group = data$population)

  q_tbl <- q_tbl %>% tidyr::pivot_longer(cols = dplyr::starts_with(".Q"),
                                         names_to = "q", values_to = "percentage")
  q_tbl$percentage <- as.numeric(q_tbl$percentage)
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
  q_tbl
}


#' Tidy ADMXITURE output files into plots
#'
#' Takes the name of a directory containing .Q file outputs, and
#' produces a list of tidied tibbles ready to plot.
#'
#' @param x the name of a directory containing .Q files
#' @param data  An associated tibble (e.g. a [`gen_tibble`]), with the individuals in the same order as the data used to
#' generate the Q matrix
#' @returns a list of `q_matrix` objects to plot
#'
#' @export

read_q_matrix_list <- function(x, data){

  files <- list.files(x, pattern = "\\.Q$", full.names = TRUE)

  # Read all .Q files into a list
  data_list <- lapply(files, function(file) utils::read.table(file, header = FALSE))

  # Turn each into a Q matrix
  matrix_list <- lapply(data_list, FUN = as_q_matrix)

  # Sort matrix_list by the number of columns
  matrix_list <- matrix_list[order(sapply(matrix_list, ncol))]

  # Tidy each
  matrix_list <- lapply(matrix_list, function(x) tidy(x, gen_tbl = data))

  matrix_list
}


#' Autoplots for `q_matrix` objects
#'
#' @param object A Q matrix object (as returned by [as_q_matrix()]).
#' @param data An associated tibble (e.g. a [`gen_tibble`]), with the individuals in the same order as the data used to
#' generate the Q matrix
#' @param annotate_group Boolean determining whether to annotate the plot with the
#' group information
#' @param ... not currently used.
#' @returns a barplot of individuals, coloured by ancestry proportion
#'
#' @export

autoplot.q_matrix <- function(object, data = NULL, annotate_group = TRUE, ...){

  rlang::check_dots_empty()
  K <- ncol(object)
  if (is.null(data)) {
    q_tbl <- as.data.frame(object)
    q_tbl$id <- 1:nrow(q_tbl)
    q_tbl <- q_tbl %>% tidyr::pivot_longer(cols = dplyr::starts_with(".Q"),
                                       names_to = "q", values_to = "percentage")
  } else {
    q_tbl <- tidy(object, data)
  }
    plt <- ggplot2::ggplot(q_tbl,
                           ggplot2::aes(x = .data$id,
                                        y = .data$percentage,
                                        fill = .data$q)) +
      ggplot2::geom_col(width = 1,
                        position = ggplot2::position_stack(reverse = TRUE))+
      ggplot2::labs(y = paste("K = ", K))+
      theme_distruct() +
      scale_fill_distruct()
    if (annotate_group){
      if (is.null(data)){
        warning("no annotation possible if 'gen_tbl' is NULL")
      } else {
        plt <- plt + annotate_group_info(q_tbl)
      }
    }
    plt
}




