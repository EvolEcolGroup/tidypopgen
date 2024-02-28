#' @export
group_by.gen_tbl <- function(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)) {
  out <- NextMethod()
  class(out) <- c("gen_grouped_df",class(out))
  return(out)

}

#' @export
dplyr_reconstruct.gen_tbl <-function(data, template)
{
  # if the genotypes are gone, drop the tbl_df class
  if (!"genotypes" %in% names(data)){
    class(data) <- class(data)[-1]
  }
  data
}
