#' Annotate an adxmiture barplot with information about the populations
#'
#' Add vertical lines to separate individuals from different populations,
#' and add optional x labels with their names
#'
#' @param q_tbl a tidied `q_matrix`
#' be in a block, all adjacent to each other
#' @returns modifier for a ggplot, added with the usual '+'
#' @keywords internal

annotate_group_info <- function(q_tbl){
  group <- q_tbl %>% dplyr::distinct(id,group) %>% dplyr::pull(group)
  if (length(rle(group)$values)!=length(unique(group))) {
    stop("values in 'group' are not ordered (they should be in consecutive blocks, one per group")
  }
# is there a way to get the aesthetic from the plot, so that we can check that
# the group variable is of the correct length???
#  if (length(group)!=object$N){
#    stop("'groups' should be of the same lenght as the original data (as found in object$N)")
#  }
  group_x <- cumsum(table(forcats::fct_inorder(group)))
  segment_data = data.frame(
    x = group_x+0.5,
    xend = group_x+0.5,
    y = rep(0,length(group_x)),
    yend = rep(1, length(group_x))
  )

  get_mid_points <- function (vec){
    (vec[-length(vec)] + vec[-1L])/2.
  }

  list(
  ggplot2::geom_segment(data = segment_data,
                        ggplot2::aes(x = .data$x,
                                     y = .data$y,
                                     xend = .data$xend,
                                     yend = .data$yend),
                        inherit.aes = FALSE),
    ggplot2::scale_x_discrete(breaks = unique(q_tbl$id)[get_mid_points(c(0,group_x))],
                                labels = names(group_x)),
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)))


}
