#' Test Hardy-Weinberg equilibrium at each locus
#'
#' Return the p-value from an exact test of HWE. It uses an implementation
#' based on [HardyWeinberg::HWExact()]), optimised for speed.
#'
#' NOTE There are no tests for this function yet! Unit tests are needed.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... further arguments to pass to [HardyWeinberg::HWExact()].
#' @returns a vector of probabilities from HWE exact test, one per locus
#' @author based on code originally written by Jan Graffleman, optimised by Andrea Manica for speed
#' @rdname loci_hwe
#' @export
loci_hwe <- function(.x, ...) {
  UseMethod("loci_hwe", .x)
}

# We should write a cpp counts function. We can't use the snp_* family of functions
# as they ignore NAs

#' @export
#' @rdname loci_hwe
loci_hwe.tbl_df <- function(.x, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_hwe(.x$genotypes, ...)
}


#' @export
#' @rdname loci_hwe
loci_hwe.vctrs_bigSNP <- function(.x, ...) {
  #rlang::check_dots_empty()
  # get the FBM
  geno_fbm <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep)>1){
    # col means for submatrix (all rows, only some columns)
    colHWE_sub <- function(X, ind, rows_to_keep) {
      apply(X[rows_to_keep, ind], 2, HWExact_fast)
    }
    hwe_p <- bigstatsr::big_apply(geno_fbm, a.FUN = colHWE_sub,
                                 rows_to_keep = rows_to_keep,
                                 ind=attr(.x,"loci")$big_index,
                                 a.combine = 'c')
  } else { # if we have a single individual
    stop ("Not implemented for a single individual")
  }
    hwe_p
}

#' @export
#' @rdname loci_hwe
loci_hwe.grouped_df <- function(.x, ...) {
  # TODO this is seriously inefficient, we need to cast it into a big_apply problem
  # of maybe it isn't that bad...
  group_map(.x, .f=~loci_hwe(.x))
}

HWExact_geno_vec <- function(x){
  HWExact_fast(c(sum(x==0, na.rm = TRUE),
    sum(x==1, na.rm = TRUE),
    sum(x==2, na.rm = TRUE)), verbose=FALSE)
}


#############
# This is a modified version of HWExact, stripped of checks to speed it up
# since we control the data that are fed to it.

HWExact_fast <- function (X, alternative = "two.sided",
                          pvaluetype = "selome",
                          eps=1e-10, x.linked = FALSE, verbose = FALSE)
{
  if(!x.linked) {
    n <- sum(X)
    Xhom <- X[c(1,3)]
    Xhet <- X[2]
    nA <- 2 * Xhom[1] + Xhet
    nB <- 2 * n - nA
    MaxHet <- min(nA, nB)
    if (MaxHet < 2) {
      pval <- 1
      prob <- 1
      pofthesample <- 1
      ind <- 1
    }
    else {
      ind <- match(Xhet, seq(MaxHet%%2, MaxHet, 2))
      enAB <- nA * nB/(2 * n - 1)
      enAB <- round(enAB, digits = 0)
      if ((enAB%%2) != (MaxHet%%2))
        enAB <- enAB + 1
      nAA <- (nA - enAB)/2
      nBB <- (nB - enAB)/2
      initialprob <- 1
      AboveExp <- NULL
      BelowExp <- NULL
      if (enAB < MaxHet)
        AboveExp <- CompProbUp(nAA, nBB, enAB, initialprob,
                                               MaxHet)
      BelowExp <- CompProbDown(nAA, nBB, enAB, initialprob)
      prob <- c(rev(BelowExp), initialprob, AboveExp)
      prob <- prob/sum(prob)
    }
    if (MaxHet%%2 == 0)
      names(prob) <- seq(0, MaxHet, 2)
    if (MaxHet%%2 == 1)
      names(prob) <- seq(1, MaxHet, 2)
    Plow <- cumsum(prob)
    Phigh <- 1 - c(0, Plow)
    Phigh <- Phigh[-length(Phigh)]
    Phwe <- pmin(1, 2 * Phigh, 2 * Plow)
    pofthesample <- prob[ind]
    pval <- switch(alternative,
                   greater = switch(pvaluetype, selome = Phigh[ind], midp = Phigh[ind] - 0.5 * pofthesample,
                                    stop("invalid value for parameter pvaluetype")),
                   less = switch(pvaluetype, selome = Plow[ind], midp = Plow[ind] - 0.5 * pofthesample,
                                 stop("invalid value for parameter pvaluetype")),
                   two.sided = switch(pvaluetype, dost = Phwe[ind], selome = sum(prob[prob <=
                                                                                        pofthesample]), midp = sum(prob[prob < pofthesample]) +
                                        0.5 * pofthesample, stop("invalid value for parameter pvaluetype")),
                   stop("invalid value for parameter alternative"))
    if (verbose) {
      D <- 0.5 * (Xhet - nA * nB/(2 * n))
      cat("Haldane Exact test for Hardy-Weinberg equilibrium (autosomal)\n")
      stringpvalue <- switch(pvaluetype, dost = "using DOST p-value\n",
                             selome = "using SELOME p-value\n", midp = "using MID p-value\n",
                             stop("invalid value for parameter pvaluetype"))
      cat(stringpvalue)
      cat(paste("sample counts: n", names(Xhom[1]), " = ",
                sep = ""), Xhom[1], paste("n", names(Xhet), " = ",
                                          sep = ""), Xhet, paste("n", names(Xhom[2]), " = ",
                                                                 sep = ""), Xhom[2], "\n")
      stringtwosided <- paste("H0: HWE (D==0), H1: D <> 0 \nD = ",
                              format(D, scientific = FALSE), "p-value = ", format(pval,
                                                                                  scientific = FALSE), "\n")
      stringgreater <- paste("H0: HWE (D==0), H1: D > 0 \nD = ",
                             format(D, scientific = FALSE), "p-value = ", format(pval,
                                                                                 scientific = FALSE), "\n")
      stringless <- paste("H0: HWE (D==0), H1: D < 0 \nD = ",
                          format(D, scientific = FALSE), "p-value = ", format(pval,
                                                                              scientific = FALSE), "\n")
      toprint <- switch(alternative, two.sided = stringtwosided,
                        greater = stringgreater, less = stringless)
      cat(toprint)
    }
  } else { # x.linked marker.
    n <- sum(X)
    nfAA <- X[1]
    nfAB <- X[2]
    nfBB <- X[3]
    nmA <- X[4]
    nmB <- X[5]
    nAf <- 2*nfAA + nfAB
    nBf <- 2*nfBB + nfAB
    nm <- nmA+nmB
    nf <- n - nm
    X <- c(nmA,nmB,nfAA,nfAB,nfBB)
    nA <- nmA + 2*nfAA + nfAB
    nB <- nmB + 2*nfBB + nfAB
    nt <- nA+nB
    pA <- nA/nt
    # immediately arrange according to minor allele
    if(nA < nB) {
      X <- c(nmA,nmB,nfAA,nfAB,nfBB)
    } else {
      X <- c(nmB,nmA,nfBB,nfAB,nfAA)
    }
    ### recompute all genotype and allele counts considering A the minor allele
    nfAA <- X[3]
    nfAB <- X[4]
    nfBB <- X[5]
    nmA <- X[1]
    nmB <- X[2]
    nAf <- 2*nfAA + nfAB
    nBf <- 2*nfBB + nfAB
    nA <- nmA + 2*nfAA + nfAB
    nB <- nmB + 2*nfBB + nfAB
    pA <- nA/nt

    Z <- auxiliartable(X)

    prob <- numeric(nrow(Z))

    pofthesample <- sample.prob.last(n,nm,nmA,nA,nfAB)

    for(i in 1:nrow(Z)) {
      prob[i] <- subsamples.prob(nA, nB, nm, nf, nt, Z[i,1], Z[i,2], Z[i,3], pofthesample, eps)
    }

    if(pvaluetype=="selome") pval <- sum(prob)
    if(pvaluetype=="midp") {
      pval <- sum(prob)-0.5*pofthesample
    }
    if(verbose) {
      cat("Graffelman-Weir exact test for Hardy-Weinberg equilibrium on the X-chromosome\n")
      stringpvalue <- switch(pvaluetype,
                             selome = "using SELOME p-value\n", midp = "using MID p-value\n",
                             stop("invalid value for parameter pvaluetype"))
      cat(stringpvalue)
      cat("Sample probability",pofthesample,"p-value = ",pval,"\n")
    }
    prob <- NA # there is no final probability vector in the x.linked case
  }
  return(pval)
}


## Unexported helper functions from HardyWeinberg
## full goes credit to the author of that package
CompProbUp <- function (nAA, nBB, EnAB, prob, MaxHet, vec = NULL)
{
  pr <- prob * 4 * nAA * nBB/((EnAB + 2) * (EnAB + 1))
  nvec <- c(vec, pr)
  if (EnAB < MaxHet - 2) {
    nvec <- CompProbUp(nAA - 1, nBB - 1, EnAB + 2, pr, MaxHet,
                       nvec)
  }
  return(nvec)
}

CompProbDown <- function (nAA, nBB, EnAB, prob, vec = NULL)
{
  pr <- prob * EnAB * (EnAB - 1)/(4 * (nAA + 1) * (nBB + 1))
  nvec <- c(vec, pr)
  if (EnAB > 3) {
    nvec <- CompProbDown(nAA + 1, nBB + 1, EnAB - 2, pr,
                         nvec)
  }
  return(nvec)
}

sample.prob.last <- function (n, nm, mA, nA, nfAB)
{
  nf <- n - nm
  nt <- nm + 2 * nf
  nB <- nt - nA
  mB <- nm - mA
  nfAA <- 0.5 * (nA - mA - nfAB)
  nfBB <- 0.5 * (nB - mB - nfAB)
  log2 <- log(2)
  K <- lgamma(nA + 1) + lgamma(nB + 1) + lgamma(nf + 1) + lgamma(nm +
                                                                   1) - lgamma(mA + 1) - lgamma(mB + 1) - lgamma(nt + 1)
  quotient <- nfAB * log2 - lgamma(nfAA + 1) - lgamma(nfAB +
                                                        1) - lgamma(nfBB + 1)
  prob <- exp(K + quotient)
  return(prob)
}

auxiliartable <- function (x, verbose = FALSE)
{
  n <- sum(x)
  nm <- x[1] + x[2]
  nf <- n - nm
  nt <- nm + 2 * nf
  nA <- x[1] + 2 * x[3] + x[4]
  nB <- nt - nA
  fAA <- x[3]
  fAB <- x[4]
  fBB <- x[5]
  mA <- x[1]
  mB <- x[2]
  if (verbose) {
    cat(n, nm, nf, "\n")
    cat(nt, nA, nB, "\n")
    cat(fAA, fAB, fBB, "\n")
    cat(mA, mB, "\n")
  }
  minor.allele.overall <- min(nA, nB)
  major.allele.overall <- max(nA, nB)
  if (verbose)
    cat("minor.allele.overall ", minor.allele.overall)
  mA.max <- min(minor.allele.overall, nm)
  mB.max <- min(major.allele.overall, nm)
  mA.min <- nm - mB.max
  mB.min <- nm - mA.max
  mA.seq <- seq(mA.min, mA.max)
  mB.seq <- seq(mB.max, mB.min, -1)
  male.minor.allele.seq.for.females <- minor.allele.overall -
    mA.seq
  male.major.allele.seq.for.females <- major.allele.overall -
    mB.seq
  minor.allele.females <- pmin(male.minor.allele.seq.for.females,
                               male.major.allele.seq.for.females)
  number.outcomes.for.females <- floor(minor.allele.females/2) +
    1
  Tab <- cbind(mA.seq, mB.seq, minor.allele.females)
  colnames(Tab) <- c("MaleMin", "MaleMaj", "FemaleMin")
  if (verbose) {
    print(Tab)
  }
  return(Tab = Tab)
}


subsamples.prob <- function (nA, nB, nm, nf, nt, mA, mB, nAf, psamp, eps)
{
  log2 <- log(2)
  nBf <- 2 * nf - nAf
  if ((nAf%%2) == 0)
    nfAB <- seq(0, nAf, 2)
  else nfAB <- seq(1, nAf, 2)
  nfAA <- 0.5 * (nAf - nfAB)
  nfBB <- 0.5 * (nBf - nfAB)
  K <- lgamma(nA + 1) + lgamma(nB + 1) + lgamma(nf + 1) + lgamma(nm +
                                                                   1) - lgamma(mA + 1) - lgamma(mB + 1) - lgamma(nt + 1)
  quotient <- nfAB * log2 - lgamma(nfAA + 1) - lgamma(nfAB +
                                                        1) - lgamma(nfBB + 1)
  prob <- exp(K + quotient)
  ii <- nearlyEqual(prob, rep(psamp, length(nfAB)), eps)
  pv <- sum(prob[ii])
  iii <- ((!ii) & (prob < psamp))
  pv <- pv + sum(prob[iii])
  return(pv)
}


nearlyEqual <- function (a, b, epsilon = 1e-10)
{
  if (length(a) != length(b))
    stop("a and b have unequal length.")
  absA <- abs(a)
  absB <- abs(b)
  diff <- abs(a - b)
  y <- rep(NA, length(a))
  i1 <- a == b
  i2 <- (a == 0 | b == 0 | (absA + absB < .Machine$double.xmin))
  y[i1] <- TRUE
  y[!i1 & i2] <- (diff[!i1 & i2] < (epsilon * .Machine$double.xmin))
  y[!i1 & !i2] <- (diff[!i1 & !i2]/min((absA[!i1 & !i2] + absB[!i1 &
                                                                 !i2]), .Machine$double.xmax) < epsilon)
  return(y)
}
