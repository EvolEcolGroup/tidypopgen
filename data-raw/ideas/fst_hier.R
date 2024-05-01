## a full break down of the fst function, which includes Nei Fst

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
test_sub_gt <- test_gt %>% filter(population%in%c("pop1","pop2"))

# convert to hierfstat
test_hier <- gt_as_hierfstat(test_sub_gt)

#hier_basic <- hierfstat::basic.stats(test_hier)
# note that Fstp, FIS and Dest are not simply averages

#hier_fst_wc <- hierfstat::pairwise.WCfst(test_hier)
#hier_fst_nei <- hierfstat::pairwise.neifst(test_hier)

########################################################
data<-test_hier
diploid = TRUE
digits=10
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

