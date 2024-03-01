#' @export
group_by.gen_tbl <- function(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)) {
  out <- NextMethod()
  class(out) <- c("grouped_df", "gen_tbl", class(out)[-1])
  return(out)

}

#' @export
dplyr_reconstruct.gen_tbl <-function(data, template)
{
  out <- NextMethod()
  # if the genotypes are gone, drop the tbl_df class
  if (!"genotypes" %in% names(data)){
    message("as genotypes were dropped, this is not longer a 'gen_tbl'")
    class(out) <- class(out)[-1]
  }
  #browser()
  out
}

#' #' @export
#' dplyr_row_slice.gen_tbl<-function(data, i, ...){
#'   #browser()
#'   NextMethod()
#'
#' }
#' #' @export
#' dplyr_col_modify.gen_tbl<-function(data, cols){
#'   #browser()
#'   NextMethod()
#'
#' }


