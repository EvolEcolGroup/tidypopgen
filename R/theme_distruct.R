#' A theme to match the output of distruct
#'
#' A theme to remove most plot decorations, matching the look of plots
#' created with distruct.
#' @returns a [ggplot2::theme]
#' @export

theme_distruct <- function(){
  ggplot2::theme_minimal()+
    ggplot2::theme( # adjust title position and remove panel grid
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank())
}
