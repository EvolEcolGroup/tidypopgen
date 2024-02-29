select_loci_if <-function(.data, .sel_logical,.col=NULL){
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  if (!inherits(.sel_bool,"logical"){
    stop(".sel_logical should be a logial (boolean) vector")
  }
  if (length(.sel_bool != ncol(show_genotypes(.data, genotypes)))){
    stop(".sel_logical should be the same length as the number of loci")
  }


}
