slidingRUNS.run <- function(x, windowSize = 15, threshold = 0.05,
                            minSNP = 3, ROHet = FALSE, maxOppWindow = 1, maxMissWindow = 1,
                            maxGap = 10^6, minLengthBps = 1000, minDensity = 1/1000,
                            maxOppRun = NULL, maxMissRun = NULL) {

  # TODO check that snps are in order (we do that for LD, I think)


  # collect all parameters in a variable
  parameters <- list(windowSize=windowSize, threshold=threshold, minSNP=minSNP,
                     ROHet=ROHet, maxOppWindow=maxOppWindow,
                     maxMissWindow=maxMissWindow, maxGap=maxGap, minLengthBps=minLengthBps,
                     minDensity=minDensity, maxOppRun=maxOppRun, maxMissRun=maxMissRun)

  # calculate gaps
  # TODO
  gaps <- diff(mapFile$bps)

  # initialize data.frame of results
  RUNs <- data.frame(group=character(), id=character(), chrom=character(), nSNP=integer(),
                     from=integer(), to=integer(), lengthBps=integer())

  chromosomes <- unique(show_loci(x)$chromosome)

  window_runs_chromosome <- function (X, ind, indiv_name, indiv_group = NULL,parameters){

  }

  # pointer to the FBM
  X <- gt_get_bigsnp(x)$genotypes

  for (i_chromosome in chromosomes){
    # get loci for this chromosome
    rows_this_chrom <- show_loci(x)$big_snp_locus




  }

  # read file line by line (http://stackoverflow.com/questions/4106764/what-is-a-good-way-to-read-line-by-line-in-r)
  while (length(oneLine <- readLines(conn, n = 1, warn = FALSE)) > 0) {

    # get individual
     individual <- list(FID=genotype[1], IID=genotype[2])

    # find runs for this individual
    a_run <- detectRUNS:::slidingRuns(genotype, animal, mapFile, gaps, parameters)

    # bind this run (if has rows) to others RUNs (if any)
    RUNs <- rbind(RUNs, a_run)

  }

  # close input stream
  close(conn)

  # fix row names
  row.names(RUNs) <- NULL

  # return calculated runs (data.frame)
  return(RUNs)
}


#' Function to detect runs using sliding window approach
#'
#' This is a core function not intended to be exported
#'
#' @param indGeno vector of 0/1/NAs of individual genotypes (0: homozygote; 1: heterozygote)
#' @param individual list of group (breed, population, case/control etc.) and ID of individual sample
#' @param mapFile Plink map file (for SNP position)
#' @param gaps distance between SNPs
#' @param parameters list of parameters
#' @param cpp use cpp functions or not (DEBUG)
#'
#' @details
#' This method uses sliding windows to detect RUNs. Checks on minimum n. of SNP, max n. of opposite and missing genotypes,
#' max gap between adjacent loci and minimum length of the run are implemented (as in the sliding window method).
#' Both runs of homozygosity (RoHom) and of heterozygosity (RoHet) can be search for (option ROHet: TRUE/FALSE)
#' NOTE: this methods is intended to not be exported
#'
#' @return A data frame of runs per individual sample
#' @keywords internal
#'

slidingRuns <- function(indGeno, individual, mapFile, gaps, parameters, cpp=TRUE) {
  # get individual and group
  ind <- as.character(individual$IID)
  group <- as.character(individual$FID)

  # use sliding windows (check cpp)
  if (cpp == TRUE) {
    res <- slidingWindowCpp(indGeno, gaps, parameters$windowSize, step=1,
                            parameters$maxGap, parameters$ROHet, parameters$maxOppWindow,
                            parameters$maxMissWindow);

    snpRun <- snpInRunCpp(res$windowStatus, parameters$windowSize, parameters$threshold)
  } else {
    res <- slidingWindow(indGeno, gaps, parameters$windowSize, step=1,
                         parameters$maxGap, parameters$ROHet, parameters$maxOppWindow,
                         parameters$maxMissWindow);

    snpRun <- snpInRun(res$windowStatus, parameters$windowSize, parameters$threshold)
  }


  # TODO: check arguments names
  dRUN <- createRUNdf(snpRun, mapFile, parameters$minSNP, parameters$minLengthBps,
                      parameters$minDensity, res$oppositeAndMissingGenotypes,
                      parameters$maxOppRun, parameters$maxMissRun)

  # manipulate dRUN to order columns
  dRUN$id <- rep(ind, nrow(dRUN))
  dRUN$group <- rep(group, nrow(dRUN))
  dRUN <- dRUN[,c(7,6,4,3,1,2,5)]

  # debug
  if(nrow(dRUN) > 0) {
    message(paste("N. of RUNS for individual", ind, "is:", nrow(dRUN)))
  } else {
    message(paste("No RUNs found for animal", ind))
  }

  #return RUNs to caller
  return(dRUN)
}

