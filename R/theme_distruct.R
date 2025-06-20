#' A theme to match the output of distruct
#'
#' A theme to remove most plot decorations, matching the look of plots
#' created with distruct.
#' @returns a [ggplot2::theme]
#' @export
#' @examples
#' # Read example gt_admix object
#' admix_obj <-
#'   readRDS(system.file("extdata", "anolis", "anole_adm_k3.rds",
#'     package = "tidypopgen"
#'   ))
#'
#' # Basic barplot with disstruct theme
#' autoplot(admix_obj, k = 3, run = 1, type = "barplot") +
#'   theme_distruct()
theme_distruct <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      # adjust title position and remove panel grid
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )
}
