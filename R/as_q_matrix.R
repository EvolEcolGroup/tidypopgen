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
    tibble::as_tibble()

  q_tbl <- q_tbl  %>%
    dplyr::mutate(id = data$id,
                  #@TODO
                  group = data$population)

  q_tbl
}


#' Augment data with information from a q_matrix object
#'
#' Augment for `q_matrix` accepts a model object and a dataset and adds
#' Q values to each observation in the dataset.
#' Q values  are stored in separate columns, which is given name with the
#' pattern ".Q1",".Q2", etc. For consistency with [broom::augment.prcomp], a column
#' ".rownames" is also returned; it is a copy of 'id', but it ensures that
#' any scripts written for data augmented with [broom::augment.prcomp] will
#' work out of the box (this is especially helpful when adapting plotting scripts).
#' @param x  A `q_matrix` object
#' @param data the `gen_tibble` used to run the clustering algorithm
#' @param ... Not used. Needed to match generic signature only.
#' @return A  [gen_tibble] containing the original data along with
#'   additional columns containing each observation's Q values.
#' @export
#' @name augment_q_matrix

augment.q_matrix <- function(x, data = NULL, ...) {

  if (!".rownames" %in% names(data)) {
    data <- data %>%
      dplyr::mutate(.rownames = data$id)
  }


  q_tbl <- tidy(x,data)

  data <- dplyr::left_join(data,q_tbl, by = "id")

}


#' Tidy ADMXITURE output files into plots
#'
#' Takes the name of a directory containing .Q file outputs, and
#' produces a list of tidied tibbles ready to plot.
#'
#' @param x the name of a directory containing .Q files
#' @returns a list of `q_matrix` objects to plot
#'
#' @export

read_q_matrix_list <- function(x){


  files <- list.files(x, pattern = "\\.Q$", full.names = TRUE)

  # Read all .Q files into a list
  data_list <- lapply(files, function(file) utils::read.table(file, header = FALSE))

  # Turn each into a Q matrix
  matrix_list <- lapply(data_list, FUN = as_q_matrix)

  # Sort matrix_list by the number of columns
  matrix_list <- matrix_list[order(sapply(matrix_list, ncol))]

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
    plt <- ggplot2::ggplot(q_tbl,
                           ggplot2::aes(x = .data$id,
                                        y = .data$percentage,
                                        fill = .data$q)) +
      ggplot2::geom_col(width = 1,
                        position = ggplot2::position_stack(reverse = TRUE))+
      ggplot2::labs(y = paste("K = ", K))+
      theme_distruct() +
      scale_fill_distruct()

    plt

    # thin vertical lines show when plot is saved as .pdf and opened with certain viewers,
    # this is a product of the specific viewer (seen on unix and mac), knitting to
    # html instead fixes, or choosing a different output (not pdf)

  } else {
    q_tbl <- tidy(object, data)

    q_tbl <- q_tbl %>% tidyr::pivot_longer(cols = dplyr::starts_with(".Q"),
                                           names_to = "q", values_to = "percentage") %>%
      dplyr::mutate(percentage = as.numeric(.data$percentage))

    q_tbl <- q_tbl %>%
      dplyr::group_by(.data$group, .data$id) %>%
      dplyr::arrange(.data$group, .data$id) %>%
      dplyr::mutate(q = factor(.data$q, levels = .data$q[order(.data$percentage, decreasing = FALSE)]))

    dominant_q <- q_tbl %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize(dominant_q = max(.data$percentage), .groups = 'drop')

    q_tbl <- q_tbl %>%
      dplyr::left_join(dominant_q, by = "id")

    q_tbl <- q_tbl %>%
      dplyr::group_by(.data$group) %>%
      dplyr::arrange(desc(.data$dominant_q), .by_group = TRUE)

    levels_q <- unique(q_tbl$id)

    q_tbl <- q_tbl %>%
      dplyr::mutate(id = factor(.data$id, levels = levels_q))

    plt <- ggplot2::ggplot(q_tbl,
                           ggplot2::aes(x = .data$id,
                                        y = .data$percentage,
                                        fill = .data$q)) +
      ggplot2::geom_col(width = 1,
                        position = "stack")+
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

}




