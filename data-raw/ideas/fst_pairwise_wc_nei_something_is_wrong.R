
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

