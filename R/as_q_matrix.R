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
#' @return A tidied matrix
#' @export
tidy.q_matrix <- function(x, gen_tbl, ...){
  rlang::check_dots_empty()
  q_tbl <- x %>%
    tibble::as_tibble() %>%
    # add the pops data for plotting
    dplyr::mutate(id = gen_tbl$id,
                  group = gen_tbl$population)

  q_tbl <- q_tbl %>% tidyr::pivot_longer(cols = dplyr::starts_with(".Q"),
                                         names_to = "q", values_to = "percentage")
  q_tbl

}
