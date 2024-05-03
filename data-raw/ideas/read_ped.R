## this is wrong. The code already reads in data as raw, so we could simply put values directly into a filebacked matrix
# if we first figure out the size (number of individuals and number of loci)


gen_tibble_ped <- function(x, ...,
                            valid_alleles = c("A", "T", "C", "G"),
                            missing_alleles = c("0","."),
                            backingfile = NULL, quiet = FALSE) {
  # Substitute .ped with .map
  map_file <- sub("\\.ped$", ".map", x)
  if (!file.exists(map_file)){
    stop("map file ",map_file," does not exist")
  }


  res <- read.pedfile(file =x , snps = map_file,
                      na.string = missing_alleles)
  # using the gen_tibble.matrix method
  new_gen_tbl <- gen_tibble(x = res$genotypes,
                            indiv_meta = res$fam,
                            loci = res$map,
                            backingfile = backingfile)
  check_allele_alphabet (new_gen_tbl, valid_alleles = valid_alleles,
                         missing_alleles = missing_alleles)
  show_loci(new_gen_tbl) <- harmonise_missing_values(show_loci(new_gen_tbl), missing_alleles = missing_alleles)
  return(new_gen_tbl)

}





# read ped using snpStats::read.pedfile

read.pedfile <- function(file, n, snps, which, split="\t| +", sep=".",
                         na.strings="0", lex.order=FALSE, quiet = FALSE) {
  ## Constants
#  r0 <- as.raw(0)
#  r1 <- as.raw(1)
#  r2 <- as.raw(2)
#  r3 <- as.raw(3)

  r0 <- 0
  r1 <- 1
  r2 <- 2
  r3 <- 3

  ## Input file
  con <- gzfile(file)
  open(con)
  ## If no line count, find out
  if (missing(n)) {
    n <- 0
    repeat {
      line <- readLines(con, n=1)
      if (length(line)==0) break;
      n  <- n+1;
    }
    if (n==0)
      stop("Nothing read")
    seek(con, 0)
  }
  ## Find snp names
  gen <- missing(snps)
  map <- NULL
  if (!gen) {
    m <- length(snps)
    if (m==1) {
      map <- read.table(snps, comment.char="")
      m <- nrow(map)
      if (missing(which)) {
        which <- 1
        repeat {
          snps <- map[,which]
          if (!any(duplicated(snps)))
            break
          if (which==ncol(map))
            stop("No unambiguous snp names found on file")
          which <- which+1
        }
      }
      else {
        snps <- map[,which]
      }
    }
  }
  else {
    line <- readLines(con, n=1)
    fields <- strsplit(line, split)[[1]]
    nf <- length(fields)
    if (nf%%2!=0)
      stop("Odd number of fields")
    m <- (nf - 6)/2
    seek(con, 0)
  }
  nf <- 6+2*m
  ## Generate empty matrix
  result <- matrix(rep(NA,n*m), nrow=n)
  ## Columns of subject dataframe
  ped <- character(n)
  mem <- character(n)
  pa <- character(n)
  ma <- character(n)
  sex <- numeric(n)
  aff <- numeric(n)
  rownms <- character(n)
  a1 <- a2 <- rep(NA, m)
  a1m <- a2m <- rep(TRUE, m)
  mallelic <- rep(FALSE, m) ## Multiallelic?
  for (i in 1:n) {
    line <- readLines(con, n=1)
    fields <- strsplit(line, "\t| +")[[1]]
    to.na <- fields %in% na.strings
    fields[to.na] <- NA
    ped[i] <- fields[1]
    mem[i] <- fields[2]
    pa[i] <- fields[3]
    ma[i] <- fields[4]
    sex[i] <- as.numeric(fields[5])
    aff[i] <- as.numeric(fields[6])
    alleles <- matrix(fields[7:nf], byrow=TRUE, ncol=2)
    one <- two <- rep(FALSE, m)
    for (k in 1:2) {
      ak <- alleles[,k]
      akm <- is.na(ak)
      br1 <- !akm & a1m
      a1[br1] <- ak[br1]
      a1m[br1] <- FALSE
      br2 <- !akm & (a1==ak)
      one[br2] <- TRUE
      br3 <- !akm & !a1m & (a1!=ak)
      br4 <- br3 & a2m
      a2[br4] <- ak[br4]
      a2m[br4] <- FALSE
      br5 <- br3 & (a2==ak)
      two[br5] <- TRUE
      mallelic <- mallelic | !(akm|one|two)
    }
    gt <- rep(r0, m)
    gt[one&!two] <- r1
    gt[one&two] <- r2
    gt[two&!one] <- r3
    result[i,] <- gt
  }
  close(con)
  ## Warnin messages
  mono <-(a2m & !a1m)

  if (any(mallelic)) {
    result[,mallelic] <- r0;
  }
  if (!quiet){
    if (any(a1m))
      message("no data for ", sum(a1m), " loci")
    if (any(mono))
      message(sum(mono), " loci were monomorphic")
    if (any(mallelic)) {
      message(sum(mallelic), " loci were multi-allelic --- set to NA")
    }
  }

  ## SnpMatrix result
  if (gen)
    snps <- paste("locus", 1:m, sep=sep)
  if (any(duplicated(ped))) {
    if (any(duplicated(mem))) {
      rnames <- paste(ped, mem, sep=sep)
      if (any(duplicated(rnames)))
        stop("could not create unique subject identifiers")
    }
    else
      rnames <- mem
  }
  else
    rnames <- ped
  dimnames(result) <- list(rnames, snps)
#  result <- new("SnpMatrix", result)
  ## Switch alleles if requested and necessary
  if (lex.order) {
    swa<- (!(is.na(a1)|is.na(a2)) & (a1>a2))
    switch.alleles(result, swa)
    a1n <- a1
    a1n[swa] <- a2[swa]
    a2[swa] <- a1[swa]
    a1 <- a1n
  }
  ## Subject support file
  fam <- data.frame(row.names=rnames, population=ped, id=mem,
                    father=pa, mother=ma, sex=sex, phenotype=aff)
  ## map data frame
  if (is.null(map))
    map <- data.frame(row.names=snps, snp.name=snps,
                      allele.1=a1, allele.2=a2)
  else {
    # mapfile
    names(map) <- c('chromosome', 'name', 'genetic_dist','position')
    map$allele_ref <- a1
    map$allele_alt <- a2
    #names(map)[which] <- "snp.names"
    map
  }
  list(genotypes=result, fam=fam, map=map)
}
