#' Annotate an admixture barplot with information about the populations
#'
#' Add vertical lines to separate individuals from different populations,
#' and add optional x labels with their names
#'
#' @param q_tbl a tidied `q_matrix`
#' be in a block, all adjacent to each other
#' @returns modifier for a ggplot, added with the usual '+'
#' @keywords internal
#' @noRd

annotate_group_info <- function(q_tbl, plt) {
  group <- q_tbl %>%
    dplyr::distinct(.data$id, .data$group) %>%
    dplyr::pull(dplyr::all_of("group"))
  if (length(rle(as.character(group))$values) != length(unique(group))) {
    stop(paste(
      "values in 'group' are not ordered",
      "(they should be in consecutive blocks, one per group"
    ))
  }
  if (length(group) != length(unique(plt$data$id))) {
    stop("'groups' should be of the same length as the original data")
  }
  group_x <- cumsum(table(fct_inorder_base(group)))
  segment_data <- data.frame(
    x = group_x + 0.5,
    xend = group_x + 0.5,
    y = rep(0, length(group_x)),
    yend = rep(1, length(group_x))
  )

  get_mid_points <- function(vec) {
    (vec[-length(vec)] + vec[-1L]) / 2.
  }

  list(
    ggplot2::geom_segment(
      data = segment_data,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        xend = .data$xend,
        yend = .data$yend
      ),
      inherit.aes = FALSE
    ),
    ggplot2::scale_x_discrete(
      breaks = unique(q_tbl$id)[get_mid_points(c(0, group_x))],
      labels = names(group_x)
    ),
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    )
  )
}
