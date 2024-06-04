#' Compute pairwise population Fst
#'
#' This function computes pairwise Fst. The following methods are implemented:
#' - 'Hudson': Hudson's formulation, as derived in Bhatia et al (2013) for diploids.
#' - 'Nei86' : Gst according to Nei (1986), as derived in Bhatia et al (2013) for diploids.
#' - 'Nei87' : Fst according to Nei (1987) - this is equivalent to `hierfstat::pairwise.neifst()`,
#' and includes the correction for heterozygosity when computing Ht
#' - 'WC84' : Weir and Cockerham (1984), as derived in Bhatia et al (2013) for diploids.
#'
#' For all formulae, the genome wide estimate is obtained by taking the ratio of the mean
#' numerators and denominators over all relevant SNPs.
#'
#' @references Bhatia G, Patterson N, Sankararaman S, Price AL. Estimating and Interpreting
#' FST: The Impact of Rare Variants. Genome Research. 2013;23(9):1514–1521.
#'
#' Nei, M. (1987) Molecular Evolutionary Genetics. Columbia University Press
#'
#' @param .x a grouped [`gen_tibble`] (as obtained by using [dplyr::group_by()])
#' @param by_locus boolean, determining whether Fst should be returned by locus(TRUE),
#' or as a single genome wide value obtained by taking the ratio of the mean numerator
#' and denominator (FALSE, the default).
#' @param method one of 'Hudson', 'Nei86', 'Nei87', and 'WC84'
#' @param n_cores number of cores to be used, it defaults to [bigstatsr::nb_cores()]
#' @returns a tibble of genome-wide pairwise Fst values with each pairwise combination
#' as a row if "by_locus=FALSE", else
#' a list including the tibble of genome-wide values as well as a matrix with pairwise
#'  Fst by locus with loci as rows and and pairwise
#' combinations as columns.
#' @export


# #' @param tidy boolean whether to return a tidy tibble. Default is TRUE, FALSE
# #' returns a matrix. THIS IS NOT IMPLEMENTED YET.


pairwise_pop_fst <- function(.x, by_locus=FALSE,
                             method= c("Hudson","Nei87","Nei86","WC84"),
                             n_cores = bigstatsr::nb_cores()){
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  method <- match.arg(method)
  if (method=="Hudson"){
    pairwise_pop_fst_hudson(.x=.x, by_locus=by_locus)
  } else if (method=="Nei87"){
    pairwise_pop_fst_nei87(.x=.x, by_locus=by_locus)
  } else if (method=="Nei86"){
    pairwise_pop_fst_nei86(.x=.x, by_locus=by_locus)
  } else if (method=="WC84"){
    pairwise_pop_fst_wc84(.x=.x, by_locus=by_locus)
  }
}

pairwise_pop_fst_hudson <- function(.x, by_locus=FALSE, n_cores = bigstatsr::nb_cores()){

  #message("this function is not properly tested yet!!!")
  #message("if you have time, test the results against something like hierfstat and create a unit test")
  message("Whilst this function should work, it has not been extensively tested. Check your results to ensure they make sense")
  # get the populations
  .group_levels = .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels),2)
  # vector and matrix to store Fst for total and by locus
  Fst_tot <- rep(NA_real_, ncol(pairwise_combn))
  if (by_locus){
    Fst_locus <- matrix(NA_real_, nrow = count_loci(.x), ncol = ncol(pairwise_combn))
  }
  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(.gt_get_bigsnp(.x)$genotypes,
                                       rowInd = .gt_bigsnp_rows(.x),
                                       colInd = .gt_bigsnp_cols(.x),
                                       groupIds = dplyr::group_indices(.x)-1,
                                       ngroups = nrow(.group_levels),
                                       ncores = n_cores)

  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]

    numerator <-  (pop_freqs_df$freq_alt[,pop1]-pop_freqs_df$freq_alt[,pop2])^2 -
      (pop_freqs_df$freq_alt[,pop1]*pop_freqs_df$freq_ref[,pop1])/ (pop_freqs_df$n[,pop1] -1) -
      (pop_freqs_df$freq_alt[,pop2]*pop_freqs_df$freq_ref[,pop2])/ (pop_freqs_df$n[,pop2] -1)
    denominator <- pop_freqs_df$freq_alt[,pop1] * pop_freqs_df$freq_ref[,pop2] +
      pop_freqs_df$freq_alt[,pop2] * pop_freqs_df$freq_ref[,pop1]
    if (by_locus){
      Fst_locus[,i_col] = numerator/denominator
    }
    Fst_tot[i_col]<-mean(numerator)/mean(denominator)
  }
  # format nicely the objects
  group_combinations <- cbind(.group_levels[pairwise_combn[1,],],.group_levels[pairwise_combn[2,],])
  names(group_combinations) <- c(paste0(dplyr::group_vars(.x),"_1"),paste0(dplyr::group_vars(.x),"_2"))
  Fst_tot <- tibble::tibble(group_combinations,value=Fst_tot)
  if (by_locus){
    rownames(Fst_locus)<-loci_names(.x)
    colnames(Fst_locus)<- apply(group_combinations,1,function(x)paste(x,collapse = "."))
  }
  if (by_locus){
    return(list(Fst_by_locus = Fst_locus, Fst = Fst_tot))
  } else{
    return(Fst_tot)
  }
}




# based on the formula in Bhatia 2013
pairwise_pop_fst_wc84 <- function(.x, by_locus=FALSE, n_cores = bigstatsr::nb_cores()){

  #message("this function is not properly tested yet!!!")
  #message("if you have time, test the results against something like hierfstat and create a unit test")
  message("Whilst this function should work, it has not been extensively tested. Check your results to ensure they make sense")
  # get the populations
  .group_levels = .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels),2)
  # vector and matrix to store Fst for total and by locus
  Fst_tot <- rep(NA_real_, ncol(pairwise_combn))
  if (by_locus){
    Fst_locus <- matrix(NA_real_, nrow = count_loci(.x), ncol = ncol(pairwise_combn))
  }
  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(.gt_get_bigsnp(.x)$genotypes,
                                       rowInd = .gt_bigsnp_rows(.x),
                                       colInd = .gt_bigsnp_cols(.x),
                                       groupIds = dplyr::group_indices(.x)-1,
                                       ngroups = nrow(.group_levels),
                                       ncores = n_cores)
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]
    n1 <- pop_freqs_df$n[,pop1]
    n2 <- pop_freqs_df$n[,pop2]
    p1 <- pop_freqs_df$freq_alt[,pop1]
    p2 <- pop_freqs_df$freq_alt[,pop2]
    q1 <- pop_freqs_df$freq_ref[,pop1]
    q2 <- pop_freqs_df$freq_ref[,pop2]
    a <- (n1*n2)/(n1+n2)
    b <- (1/(n1+n2-2))*((n1*p1*q1)+(n2*p2*q2))
    c <- (p1-p2)^2
    numerator <- 2*a*b
    denominator <- a*c+(2*a-1)*b
    if (by_locus){
      Fst_locus[,i_col] = 1-numerator/denominator
    }
    Fst_tot[i_col]<-1-mean(numerator)/mean(denominator)
  }
  # format nicely the objects
  group_combinations <- cbind(.group_levels[pairwise_combn[1,],],.group_levels[pairwise_combn[2,],])
  names(group_combinations) <- c(paste0(dplyr::group_vars(.x),"_1"),paste0(dplyr::group_vars(.x),"_2"))
  Fst_tot <- tibble::tibble(group_combinations,value=Fst_tot)
  if (by_locus){
    rownames(Fst_locus)<-loci_names(.x)
    colnames(Fst_locus)<- apply(group_combinations,1,function(x)paste(x,collapse = "."))
  }
  if (by_locus){
    return(list(Fst_by_locus = Fst_locus, Fst = Fst_tot))
  } else{
    return(Fst_tot)
  }
}

# based on the formula in Bhatia 2013
pairwise_pop_fst_nei86 <- function(.x, by_locus=FALSE, n_cores = bigstatsr::nb_cores()){

  #message("this function is not properly tested yet!!!")
  #message("if you have time, test the results against something like hierfstat and create a unit test")
  message("Whilst this function should work, it has not been extensively tested. Check your results to ensure they make sense")
  # get the populations
  .group_levels = .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels),2)
  # vector and matrix to store Fst for total and by locus
  Fst_tot <- rep(NA_real_, ncol(pairwise_combn))
  if (by_locus){
    Fst_locus <- matrix(NA_real_, nrow = count_loci(.x), ncol = ncol(pairwise_combn))
  }
  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(.gt_get_bigsnp(.x)$genotypes,
                                       rowInd = .gt_bigsnp_rows(.x),
                                       colInd = .gt_bigsnp_cols(.x),
                                       groupIds = dplyr::group_indices(.x)-1,
                                       ngroups = nrow(.group_levels),
                                       ncores = n_cores)
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]
    p1 <- pop_freqs_df$freq_alt[,pop1]
    p2 <- pop_freqs_df$freq_alt[,pop2]
    pbar <- (p1+p2)/2
    numerator <-  ((p1-p2)^2)
    denominator <-  (2*pbar*(1-pbar))
    if (by_locus){
      Fst_locus[,i_col] = numerator/denominator
    }
    Fst_tot[i_col]<-mean(numerator)/mean(denominator)
  }
  # format nicely the objects
  group_combinations <- cbind(.group_levels[pairwise_combn[1,],],.group_levels[pairwise_combn[2,],])
  names(group_combinations) <- c(paste0(dplyr::group_vars(.x),"_1"),paste0(dplyr::group_vars(.x),"_2"))
  Fst_tot <- tibble::tibble(group_combinations,value=Fst_tot)
  if (by_locus){
    rownames(Fst_locus)<-loci_names(.x)
    colnames(Fst_locus)<- apply(group_combinations,1,function(x)paste(x,collapse = "."))
  }
  if (by_locus){
    return(list(Fst_by_locus = Fst_locus, Fst = Fst_tot))
  } else{
    return(Fst_tot)
  }
}



## use tidyr::pivot_wider to turn into a matrix if that's what is requested.

### Faster versions
# the implementation for Nei 87, adapted from hierfstat
pairwise_pop_fst_nei87 <- function(.x, by_locus = FALSE, n_cores = bigstatsr::nb_cores()){
  # get the populations
  .group_levels = .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels),2)
  # vector and matrix to store Fst for total and by locus
  Fst_tot <- rep(NA_real_, ncol(pairwise_combn))
  if (by_locus){
    Fst_locus <- matrix(NA_real_, nrow = count_loci(.x), ncol = ncol(pairwise_combn))
  }
  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(.gt_get_bigsnp(.x)$genotypes,
                                       rowInd = .gt_bigsnp_rows(.x),
                                       colInd = .gt_bigsnp_cols(.x),
                                       groupIds = dplyr::group_indices(.x)-1,
                                       ngroups = nrow(.group_levels),
                                       ncores = n_cores)
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]

    n <-pop_freqs_df$n[,c(pop1,pop2)]/2
    sHo <-pop_freqs_df$het_obs[,c(pop1,pop2)]
    mHo <- apply(sHo, 1, mean, na.rm = TRUE)
    freq_alt <- pop_freqs_df$freq_alt[,c(pop1,pop2)]
    freq_ref <- pop_freqs_df$freq_ref[,c(pop1,pop2)]

    # sum of squared frequencies
    sp2 <- freq_alt^2+freq_ref^2
    Hs <- (1 - sp2 - sHo/2/n)
    Hs <- n/(n - 1) * Hs
    np <- apply(n, 1, fun <- function(x) sum(!is.na(x)))
    # mean sample size over the populations
    mn <- apply(n, 1, fun <- function(x) {
      sum(!is.na(x))/sum(1/x[!is.na(x)])
    })
    # mean sum of square frequencies
    msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
    mp2 <- rowMeans(freq_alt)^2+rowMeans(freq_ref)^2
    mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
    Ht <- 1 - mp2 + mHs/mn/np - mHo/2/mn/np

    Dst <- Ht - mHs
    Dstp <- np/(np - 1) * Dst
    Htp = mHs + Dstp
    if (by_locus){
      Fst_locus[,i_col] = Dstp/Htp
    }
    Fst_tot[i_col]<-mean(Dstp)/mean(Htp)
  }
  # format nicely the objects
  group_combinations <- cbind(.group_levels[pairwise_combn[1,],],.group_levels[pairwise_combn[2,],])
  names(group_combinations) <- c(paste0(dplyr::group_vars(.x),"_1"),paste0(dplyr::group_vars(.x),"_2"))
  Fst_tot <- tibble::tibble(group_combinations,value=Fst_tot)
  if (by_locus){
    rownames(Fst_locus)<-loci_names(.x)
    colnames(Fst_locus)<- apply(group_combinations,1,function(x)paste(x,collapse = "."))
  }
  if (by_locus){
    return(list(Fst_by_locus = Fst_locus, Fst = Fst_tot))
  } else{
    return(Fst_tot)
  }
}




