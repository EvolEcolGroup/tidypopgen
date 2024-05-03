#' Detect runs of homozygosity based on windows
#'
#' This function uses a sliding-window approach to look for runs of homozygosity (or
#' heterozygosity) in a diploid genome. This function uses the package `selectRUNS`,
#' which implements an approach equivalent to the one in PLINK.
#'
#' @param x a [gen_tibble]
#' @param window_size the size of sliding window (number of SNP loci) (default = 15)
#' @param threshold the threshold of overlapping windows of the same state
#' (homozygous/heterozygous) to call a SNP in a RUN (default = 0.05)
#' @param min_snp minimum n. of SNP in a RUN (default = 3)
#' @param heterozygosity should we look for runs of heterozygosity (instead of homozygosity? (default = FALSE)
#' @param max_opp_window max n. of SNPs of the opposite type (e.g. heterozygous
#' snps for runs of homozygosity) in the
#' sliding window (default = 1)
#' @param max_miss_window max. n. of missing SNP in the sliding window (default = 1)
#' @param max_gap max distance between consecutive SNP to be still considered a
#' potential run (default = 10^6 bps)
#' @param min_length_bps minimum length of run in bps (defaults to 1000 bps = 1 kbps)
#' @param min_density minimum n. of SNP per kbps (defaults to 0.1 = 1 SNP every 10 kbps)
#' @param max_opp_run max n. of opposite genotype SNPs in the run (optional)
#' @param max_miss_run max n. of missing SNPs in the run (optional)
#'
#' @details
#' This function returns a data frame with all runs detected in the dataset.
#' This data frame can then be written out to a csv file.
#' The data frame is, in turn, the input for other functions of the detectRUNS package
#' that create plots and produce statistics from the results
#' (see plots and statistics functions in this manual,
#' and/or refer to the detectRUNS vignette).
#'
#' If the [`gen_tibble`] is grouped, then the grouping variable is used to fill in the
#' group table. Otherwise, the group 'column' is filled with the same values as the 'id'
#' column
#'
#' @return A dataframe with RUNs of Homozygosity or Heterozygosity in the analysed dataset.
#' The returned dataframe contains the following seven columns: "group", "id", "chrom",
#' "nSNP", "from", "to", "lengthBps" (group: population, breed, case/control etc.;
#' id: individual identifier; chrom: chromosome on which the run is located;
#' nSNP: number of SNPs in the run; from: starting position of the run, in bps;
#' to: end position of the run, in bps; lengthBps: size of the run)
#' @export
#'
#' @examples
#' # don't run the example
#' if (FALSE) {
#' sheep_bed <- systems.file("extdata/sheep.bed", package="tidypopgen")
#' sheep_gt <- tidypopgen::gen_tibble(bed_path, backingfile = tempfile(), quiet=TRUE)
#' sheep_roh <- gt_roh(sheep_gt)
#' }


gt_roh <-  function(x, window_size = 15, threshold = 0.05,
                    min_snp = 3, ROHet = FALSE, max_opp_window = 1, max_miss_window = 1,
                    max_gap = 10^6, min_length_bps = 1000, min_density = 1/1000,
                    max_opp_run = NULL, max_miss_run = NULL) {
  if (!requireNamespace("detectRUNS", quietly = TRUE)) {
    stop(
      "to use this function, first install package 'detectRUNS' with\n",
      "install.packages('detectRUNS')")
  }
  # collect all parameters in a variable
  parameters <- list(windowSize=window_size, threshold=threshold, minSNP=min_snp,
                     ROHet=ROHet, maxOppWindow=max_opp_window,
                     maxMissWindow=max_miss_window, maxGap=max_gap, minLengthBps=min_length_bps,
                     minDensity=min_density, maxOppRun=max_opp_run, maxMissRun=max_miss_run)

  # create a map object
  map <- show_loci(x) %>% dplyr::select(dplyr::all_of(c("chromosome","name","position"))) %>%
    dplyr::rename("Chrom"="chromosome", "SNP"="name", "bps"="position") %>% as.data.frame()

  # compute the gaps between snps
  gaps <- diff(map$bps)
  # use groups (if defined)
  if (dplyr::is_grouped_df(sheep_gt)){
    groups <- x %>% dplyr::group_keys()
  } else {
    groups <- x$id
  }

  # initialize data.frame of results
  runs_df <- data.frame(group=character(), id=character(), chrom=character(), nSNP=integer(),
                     from=integer(), to=integer(), lengthBps=integer())
  # naively process it by row (the parallelism is implemented within individual)
  # access time is horrible, but I don't think this is the bottleneck
  # it needs some profiling
  X <- .gt_get_bigsnp(x)$genotypes # pointer for the FBM
  for (i in seq_len(nrow(x))){
    this_genotype <- X[i,]
    this_indiv <- list(FID=groups[i], IID=x$id[i])
    # find runs for this individual
    this_runs <- utils::getFromNamespace("slidingRuns", "detectRUNS")(this_genotype, this_indiv, map, gaps, parameters)
    # bind this run (if has rows) to others RUNs (if any)
    runs_df <- rbind(runs_df, this_runs)
  }

  # fix row names
  row.names(runs_df) <- NULL

  # return calculated runs (data.frame)
  return(runs_df)

}
