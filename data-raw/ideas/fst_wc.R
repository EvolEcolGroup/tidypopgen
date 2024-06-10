test_genotypes <- rbind(c(1,1,0,1,1,0),
                        c(2,1,0,NA,0,0),
                        c(2,NA,0,0,1,1),
                        c(1,0,0,1,0,0),
                        c(1,2,0,1,2,1),
                        c(0,0,0,0,NA,1),
                        c(0,1,1,0,1,NA))
test_indiv_meta <- data.frame (id=c("a","b","c","d","e","f","g"),
                               population = c("pop1","pop1","pop2","pop2","pop1","pop3","pop3"))
test_loci <- data.frame(name=paste0("rs",1:6),
                        chromosome=paste0("chr",c(1,1,1,1,2,2)),
                        position=as.integer(c(3,5,65,343,23,456)),
                        genetic_dist = as.integer(rep(0,6)),
                        allele_ref = c("A","T","C","G","C","T"),
                        allele_alt = c("T","C", NA,"C","G","A"))

test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)
test_gt <- test_gt %>% group_by(population)

# convert to hierfstat
test_hier <- gt_as_hierfstat(test_gt)

hier_basic <- hierfstat::basic.stats(test_hier)
# note that Fstp, FIS and Dest are not simply averages

hier_fst_wc <- hierfstat::pairwise.WCfst(test_hier)
hier_fst_nei <- hierfstat::pairwise.neifst(test_hier)

# a developer function to create various count summaries of a population, used to
# compute more complex statistics (e.g. pairwise fst, etc.).
# Specifically, we compute:
# sums_alt sum alleles for alt
# sums_ref sums alleles for ref
# n sums of all alleles
# het_obs observed heterozygosity
# het_exp expected heterozygosity
.gt_pop_counts <- function(.x){
  counts <- bigstatsr::big_counts( .gt_get_bigsnp(.x)$genotypes,
                                   ind.row =.gt_bigsnp_rows(.x),
                                   ind.col = .gt_bigsnp_cols(.x))
  sums_alt <- apply(counts,2,function(x) x[2]+2*x[3])
  n <- apply(counts,2,function(x) sum(x[1:3])*2)
  sums_ref <- n - sums_alt
  freq_alt <- sums_alt/n
  freq_ref <- 1- freq_alt
  het_count <- apply(counts,2,function(x) x[2]) # This should be simply a row from the counts matrix
  het_obs <- apply(counts,2,function(x) x[2]/sum(x[1:3]))
  het_exp <- 2 * sums_alt/n * sums_ref/n
  return (list(sums_alt = sums_alt,
               sums_ref = sums_ref,
               freq_alt = freq_alt,
               freq_ref = freq_ref,
               het_count = het_count,
               n = n,
               het_obs = het_obs,
               het_exp = het_exp))
}

pairwise_pop_fst_nei73 <- function(.x, by_locus = FALSE){
  pop_counts_df <- group_map(.x, .f=~.gt_pop_counts(.x))
  # get the grouping column, and create all pairwise combination of indices
  .group_levels = .x %>% group_keys()
  pairwise_combn <- utils::combn(nrow(.group_levels),2)
  n_loci <- count_loci(.x)
  Hs_pair <- matrix(NA_real_, nrow = n_loci, ncol = ncol(pairwise_combn))
  Ht_pair <- matrix(NA_real_, nrow = n_loci, ncol = ncol(pairwise_combn))
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]
    Hs_pair[,i_col] <-  (pop_counts_df[[pop1]]$het_exp + pop_counts_df[[pop2]]$het_exp)/2
    # Ht_pair is 1-p_bar^2 - q_bar^2
    Ht_pair[, i_col] <-
      1 - ((pop_counts_df[[pop1]]$sums_alt + pop_counts_df[[pop2]]$sums_alt) /
             (pop_counts_df[[pop1]]$n + pop_counts_df[[pop2]]$n)) ^ 2 -
      ((pop_counts_df[[pop1]]$sums_ref + pop_counts_df[[pop2]]$sums_ref) /
         (pop_counts_df[[pop1]]$n + pop_counts_df[[pop2]]$n)) ^ 2
  }
  if (by_locus){
    1-(Hs_pair/Ht_pair)
    # tidy it properly
    #TODO
  } else {
    1-(colSums(Hs_pair,na.rm=TRUE)/colSums(Ht_pair,na.rm=TRUE))
    # TODO TIDY IT PROPERLY
  }
}

pairwise_pop_fst_nei83 <- function(.x, by_locus = FALSE){
  pop_counts_df <- group_map(.x, .f=~.gt_pop_counts(.x))
  # get the grouping column, and create all pairwise combination of indices
  .group_levels = .x %>% group_keys()
  pairwise_combn <- utils::combn(nrow(.group_levels),2)
  n_loci <- count_loci(.x)
#  Hs_pair <- matrix(NA_real_, nrow = n_loci, ncol = ncol(pairwise_combn))
#  Ht_pair <- matrix(NA_real_, nrow = n_loci, ncol = ncol(pairwise_combn))
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]
    Hs_pair <-  (pop_counts_df[[pop1]]$het_count + pop_counts_df[[pop2]]$het_count)/
      (pop_counts_df[[pop1]]$n + pop_counts_df[[pop2]]$n)
    # Ht_pair is 1-p_bar^2 - q_bar^2
    Ht_pair[, i_col] <-
      1 - ((pop_counts_df[[pop1]]$sums_alt + pop_counts_df[[pop2]]$sums_alt) /
             (pop_counts_df[[pop1]]$n + pop_counts_df[[pop2]]$n)) ^ 2 -
      ((pop_counts_df[[pop1]]$sums_ref + pop_counts_df[[pop2]]$sums_ref) /
         (pop_counts_df[[pop1]]$n + pop_counts_df[[pop2]]$n)) ^ 2
    Nar <- n <- 2.0 / ((1.0 / (pop_counts_df[[pop1]]$n/2.0)) + (1.0 / (pop_counts_df[[pop2]]$n/2.0)))
    Hs_pair_est=Hs_pair*((2.0*Nar)/(2.0*Nar-1))
    Ht_pair_est=Ht_pair+(Hs_pair_est/(4.0*Nar))
  }
  if (by_locus){
    1-(Hs_pair_est/Ht_pair_est)
    # tidy it properly
    #TODO
  } else {
    1-(colSums(Hs_pair_est,na.rm=TRUE)/colSums(Ht_pair_est,na.rm=TRUE))
    # TODO TIDY IT PROPERLY
  }
}



pairwise_pop_fst_nei73(test_gt)

pop1<- test_gt %>%filter(population=="pop1")
pop1_hier <- gt_as_hierfstat(pop1)

test_gt%>%group_by(population)->.x
pop_counts_df <- group_map(.x, .f=~.gt_pop_counts(.x))

hierfstat::basic.stats(test_hier)->foo

########################################################
data<-test_sub_hier
diploid = TRUE
dum.pop<-FALSE
if (length(table(data[, 1])) < 2){
  data[dim(data)[1] + 1, 1] <- "DumPop"
  dum.pop<-TRUE
}
if (dim(data)[2] == 2)
  data <- data.frame(data, dummy.loc = data[, 2])
loc.names <- names(data)[-1]
p <- hierfstat:::pop.freq(data, diploid)
n <- t(hierfstat:::ind.count(data))
if (diploid) {
  dum <- hierfstat::getal.b(data[, -1])
  Ho <- dum[, , 1] == dum[, , 2]
  sHo <- (1 - t(apply(Ho, 2, fun <- function(x) tapply(x,
                                                       data[, 1], mean, na.rm = TRUE))))
  mHo <- apply(sHo, 1, mean, na.rm = TRUE)
} else {
  sHo <- NA
  mHo <- NA
}
sp2 <- lapply(p, fun <- function(x) apply(x, 2, fun2 <- function(x) sum(x^2))) # cols sums of squares
sp2 <- matrix(unlist(sp2), nrow = dim(data[, -1])[2], byrow = TRUE)
if (diploid) {
  Hs <- (1 - sp2 - sHo/2/n)
  Hs <- n/(n - 1) * Hs
  Fis = 1 - sHo/Hs
} else {
  Hs <- n/(n - 1) * (1 - sp2)
  Fis <- NA
}
np <- apply(n, 1, fun <- function(x) sum(!is.na(x)))
mn <- apply(n, 1, fun <- function(x) {
  np <- sum(!is.na(x))
  np/sum(1/x[!is.na(x)])
})
msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
mp <- lapply(p, fun <- function(x) apply(x, 1, mean, na.rm = TRUE))
mp2 <- unlist(lapply(mp, fun1 <- function(x) sum(x^2)))
if (diploid) {
  mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
  Ht <- 1 - mp2 + mHs/mn/np - mHo/2/mn/np
  mFis = 1 - mHo/mHs
} else {
  mHs <- mn/(mn - 1) * (1 - msp2)
  Ht <- 1 - mp2 + mHs/mn/np
  mFis <- NA
}
Dst <- Ht - mHs
Dstp <- np/(np - 1) * Dst
Htp = mHs + Dstp
Fst = Dst/Ht
Fstp = Dstp/Htp
Dest <- Dstp/(1 - mHs)
res <- data.frame(cbind(mHo, mHs, Ht, Dst, Htp, Dstp, Fst,
                        Fstp, mFis, Dest))
names(res) <- c("Ho", "Hs", "Ht", "Dst", "Htp", "Dstp",
                "Fst", "Fstp", "Fis", "Dest")
if (diploid) {
  rownames(sHo) <- loc.names
  rownames(Fis) <- loc.names
}
is.na(res) <- do.call(cbind, lapply(res, is.infinite))
overall <- apply(res, 2, mean, na.rm = TRUE)
overall[7] <- overall[4]/overall[3]
overall[8] <- overall[6]/overall[5]
overall[9] <- 1 - overall[1]/overall[2]
overall[10] <- overall[6]/(1 - overall[2])
names(overall) <- names(res)
if (!diploid) {
  overall[c(1, 9)] <- NA
}
if(dum.pop){
  ToBeRemoved<-which(dimnames(Hs)[[2]]=="DumPop")
  n<-n[,-ToBeRemoved,drop=FALSE]
  Hs<-Hs[,-ToBeRemoved,drop=FALSE]
  if (diploid) sHo<-sHo[,-ToBeRemoved,drop=FALSE]
  Fis<-Fis[,-ToBeRemoved,drop=FALSE]
  p<-lapply(p,function(x) x[,-ToBeRemoved,drop=FALSE])
}
all.res <- list(n.ind.samp = n, pop.freq = lapply(p, round,
                                                  digits), Ho = round(sHo, digits), Hs = round(Hs, digits),
                Fis = round(Fis, digits), perloc = round(res, digits),
                overall = round(overall, digits))
class(all.res) <- "basic.stats"
all.res



## steps
# create a hierfstat object with just two populations
test_sub_gt <- test_gt %>% filter(population%in%c("pop1","pop2"))
# contrast our estimates with Jerome's estimates, making sure that we use all the appropriate corrections as implemented in hierfstat
test_sub_hier <- test_sub_gt %>% gt_as_hierfstat()
# rename the dataset to the internal name in the functions to follow the estimates
data <- test_sub_hier
.x <- test_sub_gt %>% group_by(population)

n <- do.call(cbind, lapply(pop_counts_df, function(X) X[["n"]]))
sHo <- do.call(cbind, lapply(pop_counts_df, function(X) X[["het_obs"]]))
# sum of squared frequencies
sp2 <- lapply(pop_counts_df, fun <- function(x) ((x$freq_alt)^2+(x$freq_ref)^2))
sp2 <- matrix(unlist(sp2, use.names = FALSE), ncol=2, byrow = FALSE)

# square of mean frequency
mp <- lapply(p, fun <- function(x) apply(x, 1, mean, na.rm = TRUE))
sp2 <- lapply(pop_counts_df, fun <- function(x) ((x$freq_alt)+(x$freq_ref))/2)

mp2 <- unlist(lapply(mp, fun1 <- function(x) sum(x^2)))

Hs <- (1 - sp2 - sHo/2/n)
Hs <- n/(n - 1) * Hs

np <- apply(n, 1, fun <- function(x) sum(!is.na(x)))
mn <- apply(n, 1, fun <- function(x) {
  np <- sum(!is.na(x))
  np/sum(1/x[!is.na(x)])
})
msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
mp <- lapply(p, fun <- function(x) apply(x, 1, mean, na.rm = TRUE))
mp2 <- unlist(lapply(mp, fun1 <- function(x) sum(x^2)))
mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
Ht <- 1 - mp2 + mHs/mn/np - mHo/2/mn/np

Dst <- Ht - mHs
Dstp <- np/(np - 1) * Dst
Htp = mHs + Dstp
Fst = Dst/Ht
