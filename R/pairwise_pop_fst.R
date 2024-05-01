#' Compute pairwise population Fst
#'
#' This function computes pairwise Fst. The following methods are implemented:
#' - 'Hudson': Hudson's formulation, as derived in Bhatia et al (2013) for diploids.
#' - 'Nei86' : Gst according to Nei (1986), as derived in Bhatia et al (2013) for diploids.
#' - 'Nei87': Fst according to Nei (1987) - this is equivalent to `hierfstat::pairwise.neifst()`
#' - 'WC84' : Weir adn Cockerham (1984), as derived in Bhatia et al (2013) for diploids.
#'
#' For all formulae, the genome wide estimate is obtained by taking the ratio of the mean
#' numerators and denominators over all relevant SNPs.
#'
#' @references Bhatia G, Patterson N, Sankararaman S, Price AL. Estimating and Interpreting
#' FST: The Impact of Rare Variants. Genome Research. 2013;23(9):1514–1521.
#'
#' Nei, M. (1987) Molecular Evolutionary Genetics. Columbia University Press
#' @param .x a grouped [`gen_tibble`] (as obtained by using [dplyr::group_by()])
#' @param by_locus boolean, determining whether Fst should be returned by locus(TRUE),
#' or as a single genome wide value obtained by taking the ratio of the mean numerator
#' and denominator (FALSE, the default).
#' @param method one of 'Hudson', 'Nei86', 'Nei87', and 'WC84'
#' @returns a tibble of genome-wide pairwise Fst values with each pairwise combination
#' as a row if "by_locus=FALSE", else
#' a list including the tibble of genome-wide values as well as a matrix with pairwise
#'  Fst by locus with loci as rows and and pairwise
#' combinations as columns.
#' a matrix
#' @export


# #' @param tidy boolean whether to return a tidy tibble. Default is TRUE, FALSE
# #' returns a matrix. THIS IS NOT IMPLEMENT YET.


pairwise_pop_fst <- function(.x, by_locus=FALSE, method= c("Hudson","Nei87")){
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  method <- match.arg(method)
  if (method=="Hudson"){
    pairwise_pop_fst_hudson(.x=.x, by_locus=by_locus)
  } else if (method=="Nei87"){
    pairwise_pop_fst_nei87(.x=.x, by_locus=by_locus)
  }
}

pairwise_pop_fst_hudson <- function(.x, by_locus=FALSE){

  message("this function is not properly tested yet!!!")
  message("if you have time, test the results against something like hierfstat and create a unit test")
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
  pop_freqs_df <- group_map(.x, .f=~.gt_pop_freqs(.x))
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]

    numerator <-  (pop_freqs_df[[pop1]]$freq_alt-pop_freqs_df[[pop2]]$freq_alt)^2 -
      (pop_freqs_df[[pop1]]$freq_alt*pop_freqs_df[[pop1]]$freq_ref)/ (pop_freqs_df[[pop1]]$n -1) -
      (pop_freqs_df[[pop2]]$freq_alt*pop_freqs_df[[pop2]]$freq_ref)/ (pop_freqs_df[[pop2]]$n -1)
    denominator <- pop_freqs_df[[pop1]]$freq_alt * pop_freqs_df[[pop2]]$freq_ref +
      pop_freqs_df[[pop2]]$freq_alt * pop_freqs_df[[pop1]]$freq_ref
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
    rownames(Fst_locus)<-show_loci_names(.x)
    colnames(Fst_locus)<- apply(group_combinations,1,function(x)paste(x,collapse = "."))
  }
  if (by_locus){
    return(list(Fst_by_locus = Fst_locus, Fst = Fst_tot))
  } else{
    return(Fst_tot)
  }
}




# the implementation for Nei 87, adapted from hierfstat
pairwise_pop_fst_nei87 <- function(.x, by_locus = FALSE){
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
  pop_freqs_df <- group_map(.x, .f=~.gt_pop_freqs(.x))
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]

    n <-cbind(pop_freqs_df[[pop1]]$n,pop_freqs_df[[pop2]]$n)/2
    sHo <-cbind(pop_freqs_df[[pop1]]$het_obs,pop_freqs_df[[pop2]]$het_obs)
    mHo <- apply(sHo, 1, mean, na.rm = TRUE)
    freq_alt <- cbind(pop_freqs_df[[pop1]]$freq_alt, pop_freqs_df[[pop2]]$freq_alt)
    freq_ref <- cbind(pop_freqs_df[[pop1]]$freq_ref, pop_freqs_df[[pop2]]$freq_ref)

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
    rownames(Fst_locus)<-show_loci_names(.x)
    colnames(Fst_locus)<- apply(group_combinations,1,function(x)paste(x,collapse = "."))
  }
  if (by_locus){
    return(list(Fst_by_locus = Fst_locus, Fst = Fst_tot))
  } else{
    return(Fst_tot)
  }
}


# based on the formula in Bhatia 2013
pairwise_pop_fst_wc <- function(.x, by_locus=FALSE){

  message("this function is not properly tested yet!!!")
  message("if you have time, test the results against something like hierfstat and create a unit test")
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
  pop_freqs_df <- group_map(.x, .f=~.gt_pop_freqs(.x))
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]
    n1 <- pop_freqs_df[[pop1]]$n
    n2 <- pop_freqs_df[[pop2]]$n
    p1 <- pop_freqs_df[[pop1]]$freq_alt
    p2 <- pop_freqs_df[[pop2]]$freq_alt
    q1 <- pop_freqs_df[[pop1]]$freq_ref
    q2 <- pop_freqs_df[[pop2]]$freq_ref
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
    rownames(Fst_locus)<-show_loci_names(.x)
    colnames(Fst_locus)<- apply(group_combinations,1,function(x)paste(x,collapse = "."))
  }
  if (by_locus){
    return(list(Fst_by_locus = Fst_locus, Fst = Fst_tot))
  } else{
    return(Fst_tot)
  }
}

# based on the formula in Bhatia 2013
pairwise_pop_fst_nei86 <- function(.x, by_locus=FALSE){

  message("this function is not properly tested yet!!!")
  message("if you have time, test the results against something like hierfstat and create a unit test")
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
  pop_freqs_df <- group_map(.x, .f=~.gt_pop_freqs(.x))
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]
    p1 <- pop_freqs_df[[pop1]]$freq_alt
    p2 <- pop_freqs_df[[pop2]]$freq_alt
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
    rownames(Fst_locus)<-show_loci_names(.x)
    colnames(Fst_locus)<- apply(group_combinations,1,function(x)paste(x,collapse = "."))
  }
  if (by_locus){
    return(list(Fst_by_locus = Fst_locus, Fst = Fst_tot))
  } else{
    return(Fst_tot)
  }
}



## use tidyr::pivot_wider to turn into a matrix if that's what is requested.




