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
  .data$genotypes <- lapply(.data$genotypes, .SNPbin_subset, .sel_logical)
  .data
}


.SNPbin_subset <- function(x, i){
  if (missing(i)) i <- TRUE
  temp <- .SNPbin2int(x) # data as integers with NAs
  x <- new("SNPbin", snp=temp[i], label=x@label, ploidy=x@ploidy)
  return(x)
}

#############
## .SNPbin2int
#############
## convert SNPbin to integers (0/1/2...)
.SNPbin2int <- function(x){
  ##res <- lapply(x@snp, .raw2bin)
  resSize <- length(x@snp[[1]])*8
  # Wed Apr 12 08:49:02 2017 ------------------------------
  # I am leaving this function along as it does not necessarily break solaris,
  # but I am leaving the code and timings just in case.
  #
  # ZNK
  res <- .C("bytesToInt", unlist(x@snp), length(x@snp[[1]]), length(x@snp),
            integer(resSize), as.integer(resSize), PACKAGE="adegenet")[[4]][1:nLoc(x)]
  # library(microbenchmark)
  # set.seed(5000)
  # dat <- sample(c(0:2,NA), 1e5, prob=c(rep(.995/5,3), 0.005), replace=TRUE)
  # x <- new("SNPbin", dat)
  # y <- microbenchmark(C = .SNPbin2int(x), base = .SNPbin2int1(x), times = 1000)
  # print(y, "relative")
  ## Unit: relative
  ##  expr      min       lq     mean   median      uq       max neval cld
  ##     C 1.000000 1.000000 1.000000 1.000000 1.00000 1.0000000  1000   a
  ##  base 2.206831 1.200831 1.019783 1.168149 1.02654 0.1668163  1000   a
  # res     <- vapply(x@snp, function(x) as.integer(rawToBits(x)), integer(resSize))
  # res     <- apply(res[1:nLoc(x), ], 1, sum)
  if (length(x@NA.posi) > 0){
    res[x@NA.posi] <- NA_integer_
  }
  return(res)
} # end .SNPbin2int
