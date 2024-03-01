#' The `select_if` verb for `loci`
#'
#' An equivalent to [dplyr::select_if()] that works on the `genotype` column
#' of a `gen_tibble`. This function has access to the genotypes (and thus can
#' work on summary statistics to select), but not the names of the loci (look
#' at [select_loci()] for that feature. There is also the possibility
#' of swapping reference and alternate alleles whilst selecting loci, using
#' either the argument `.swap_arg` (which has `select`-like mini-grammar on names)
#' or `.swap_arg_if` (which has access to the genotypes like `select_loci_if`)
#' @param .data a `gen_tibble`
#' @param .sel_logical a logical vector of length equal to the number of loci,
#' or an expression that will tidy evaluate to such a vector.
#' @param .swap_arg one unquoted expression, using the mini-grammar of
#'  [dplyr::select()] determining at which loci the reference and alterante
#'  allele will be swapped.
#' @param .swap_if_arg a logical vector of length equal to the number of loci,
#' or an expression that will tidy evaluate to such a vector.
#' @returns a list of `SNPbin` object that have been sub-setted to the selected
#' loci.
#' @export
#'
select_loci_if <-function(.data, .sel_logical, .swap_arg = NULL, .swap_if_arg = NULL){
  # defuse the boolean argument
  sel_defused <- rlang::enquo(.sel_logical)
  # defuse the swap criteria
  swap_defused <- rlang::enquo(.swap_arg)
  swap_if_defused <- rlang::enquo(.swap_if_arg)

  if (all(!rlang::quo_is_null(swap_defused),!rlang::quo_is_null(swap_if_defused))){
    stop("only one of 'swap' or 'swap_if' can be used at any one time")
  }

  # and now evaluate it, allowing it to see the data
  loci_sel <- rlang::eval_tidy(sel_defused,data=.data)
  if (!inherits(loci_sel,"logical")){
    stop(".sel_logical should be a logical (boolean) vector")
  }
  if (length(loci_sel) != ncol(show_genotypes(.data$genotypes))){
    stop(".sel_logical should be the same length as the number of loci")
  }
  #extract the loci table
  loci_info <- attr(.data$genotypes,"loci")

  # subset the genotypes
  if (all(rlang::quo_is_null(swap_defused),rlang::quo_is_null(swap_if_defused))){
    .data$genotypes <- lapply(.data$genotypes, .SNPbin_subset, loci_sel)
  } else{ # subset and swap genotypes (slower)

    # if we also want to swap based on a select argument
    if (!rlang::quo_is_null(swap_defused)){
      # create a named vector of all loci
      loci_names_vec <- stats::setNames(seq_len(nrow(show_loci(.data))),
                                        nm = show_loci_names(.data))
      loci_swap <- tidyselect::eval_select(expr=swap_defused, data = loci_names_vec)
    } else { # based on an if (boolean) argument
      loci_swap <- rlang::eval_tidy(swap_if_defused,data=.data)
    }
    # swap alleles in the snpbin_list
    .data$genotypes <- lapply(.data$genotypes, select_swap_snpbin,
                              sel_indices = loci_sel,
                              to_swap = loci_swap)
    # swap alleles in table
    old_ref <- loci_info[loci_swap,"allele_ref"]
    loci_info[loci_swap,"allele_ref"] <- loci_info[loci_swap,"allele_alt"]
    loci_info[loci_swap,"allele_alt"] <- old_ref
  }
  # reintroduce the subsetted table
  attr(.data$genotypes,"loci") <- loci_info[loci_sel,]
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
