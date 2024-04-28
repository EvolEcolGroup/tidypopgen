#' Compute pairwise population Fst
#'
#' This function computes pairwise Fst, currently using the formulation by
#' Hudson (following the recommendations by Bhatia et al (2013)).
#'
#' @references Bhatia G, Patterson N, Sankararaman S, Price AL. Estimating and Interpreting
#' FST: The Impact of Rare Variants. Genome Research. 2013;23(9):1514–1521.
#' @param .x a grouped [`gen_tibble`] (as obtained by using [dplyr::group_by()])
#' @param by_locus boolean, determining whether Fst should be returned by locus(TRUE),
#' or as a single genome wide value obtained by taking the ratio of the mean numerator
#' and denominator (FALSE, the default).
#' @param method currently only "Hudson" is implemented
#' @returns a tibble with each pairwise combination as a row (if `tidy=TRUE`), else
#' a matrix
#' @export


# #' @param tidy boolean whether to return a tidy tibble. Default is TRUE, FALSE
# #' returns a matrix. THIS IS NOT IMPLEMENT YET.


pairwise_pop_fst <- function(.x, by_locus=FALSE, method= "Hudson"){
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  method <- match.arg(method)
  if (method=="Hudson"){
    pairwise_pop_fst_hudson(.x=.x, by_locus=by_locus)
  }
}


pairwise_pop_fst_hudson <- function(.x, by_locus=FALSE){

  message("this function is not properly tested yet!!!")
  message("if you have time, test the results against something like hierfstat and create a unit test")
  # matrix(unlist(z, use.names = FALSE), ncol = 10, byrow = TRUE)
  # is known to be faster than do.call(rbind,f)
  # see https://stackoverflow.com/questions/13224553/how-to-convert-a-huge-list-of-vector-to-a-matrix-more-efficiently
  # do not modify the approach used below:
  n_loci <- nrow(show_loci(.x))
  # sum alt alleles over each locus in each group
  sums <- matrix(unlist(.x %>%
                          loci_sums(),
                        use.names = FALSE), ncol = n_loci, byrow = TRUE)
  # get the total  number of alleles (i.e. removing NAs) for each locus in each group
  n_missing <- matrix(unlist(.x %>%
                       loci_missingness(as_counts=TRUE),
                     use.names = FALSE), ncol = n_loci, byrow = TRUE)
  n_ind <- (.x %>% tally() %>% pull(.data$n)) # number of individuals per population

  n_alleles <- sweep(-n_missing,1,n_ind, FUN = "+")*2
  # function to compute het by row
  het_exp_by_row <- function(i, sums, n){(sums[i,]*(n[i,]-sums[i,]))/(n[i,]*(n[i,]-1))}
  # get het at each locus for each population
  het <- matrix(unlist(lapply(1:nrow(sums), het_exp_by_row, sums, n_alleles),
                       use.names = FALSE), ncol = n_loci, byrow = TRUE)
  # get the grouping column, and create all pairwise combination of indices
  .group_levels = .x %>% group_keys()
  pairwise_combn <- t(utils::combn(nrow(.group_levels),2))
  numerator <- matrix(NA_real_, nrow = nrow(pairwise_combn), ncol = n_loci)
  denominator <- matrix(NA_real_, nrow = nrow(pairwise_combn), ncol = n_loci)
  # TODO this could be made faster by parallelising it
  for (i_row in seq_len(nrow(pairwise_combn))){
    pop1 <- pairwise_combn[i_row,1]
    pop2 <- pairwise_combn[i_row,2]
    numerator[i_row,] <- (sums[pop1,]/n_alleles[pop1,] - sums[pop2,]/n_alleles[pop2,])^2 -
      het[pop1,]/n_alleles[pop1,] - het[pop2,]/n_alleles[pop2,]
    denominator[i_row,] <- numerator[i_row,] + het[pop1,] + het[pop2,]
  }

  # set up first two columns with pairs of populations
  fst <- tibble(!!paste0(names(.group_levels),"_",1) := .group_levels[pairwise_combn[,1],1] %>% pull(1),
                !!(paste0(names(.group_levels),"_",2)) := .group_levels[pairwise_combn[,2],1] %>% pull(1))
  # if we want whole genome level estimates
  if (!by_locus){
    fst_wg_by_row <- function(i, numerator, denominator){
      return(sum(numerator[i,],na.rm=TRUE)/sum(denominator[i,],na.rm=TRUE))
    }
    fst <- fst %>% mutate(value = unlist(lapply(seq_len(nrow(fst)),fst_wg_by_row, numerator, denominator)))
    return(fst)
  } else {
    # if we want estimates by locus
    fst_by_row <- function(i, numerator, denominator){
      return(numerator[i,]/denominator[i,])
    }
    fst_vals <- matrix(unlist(lapply(seq_len(nrow(numerator)),fst_by_row, numerator, denominator),
                              use.names = FALSE), ncol = n_loci, byrow = TRUE)
    colnames(fst_vals) <- show_loci_names(.x)
    return(fst %>% cbind(fst_vals))
  }
}

## use tidyr::pivot_wider to turn into a matrix if that's what is requested.

