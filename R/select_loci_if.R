#' The `select_if` verb for `loci`
#'
#' An equivalent to [dplyr::select_if()] that works on the `genotype` column
#' of a `gen_tibble`. Tidy evaluation should work as expected.
#' @param .data a `gen_tibble`
#' @param .sel_logical a logical vector of length equal to the number of loci,
#' or an expression that will tidy evaluate to such a vector
#' @returns a list of `SNPbin` object that have been subsetted.
#' @export
#'
select_loci_if <-function(.data, .sel_logical){
  # defuse the boolean argument
  sel_defused <- rlang::enquo(.sel_logical)
  # and now evaluate it, allowing it to see the data
  .sel_logical <- rlang::eval_tidy(sel_defused,data=.data)
  if (!inherits(.sel_logical,"logical")){
    stop(".sel_logical should be a logical (boolean) vector")
  }
  if (length(.sel_logical) != ncol(show_genotypes(.data$genotypes))){
    stop(".sel_logical should be the same length as the number of loci")
  }
  #TODO we need to get the loci table and move it over
  loci_info <- attr(.data$genotypes,"loci")
  .data$genotypes <- lapply(.data$genotypes, .SNPbin_subset, .sel_logical)
  attr(.data$genotypes,"loci") <- loci_info[.sel_logical,]
  .data
}


# this a copy of an unexported function from adegenet
.SNPbin_subset <- function(x, i){
  if (missing(i)) i <- TRUE
  temp <- .SNPbin2int(x) # data as integers with NAs
  x <- methods::new("SNPbin", snp=temp[i], label=x@label, ploidy=x@ploidy)
  return(x)
}

#############
## .SNPbin2int
#############
## convert SNPbin to integers (0/1/2...)
# this a copy of an unexported function from adegenet
.SNPbin2int <- function(x){
  resSize <- length(x@snp[[1]])*8
  res <- .C("bytesToInt", unlist(x@snp), length(x@snp[[1]]), length(x@snp),
            integer(resSize), as.integer(resSize), PACKAGE="adegenet")[[4]][1:nLoc(x)]
  if (length(x@NA.posi) > 0){
    res[x@NA.posi] <- NA_integer_
  }
  return(res)
}
