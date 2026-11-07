#' Detect runs of homozygosity using a sliding-window approach
#'
#' This function uses a sliding-window approach to look for runs of homozygosity
#' (or heterozygosity) in a diploid genome. It implements the algorithm of
#' [detectRUNS::slidingRUNS.run()], which is equivalent to the one in PLINK,
#' natively in C++: all individuals are analysed in a single pass over the
#' genotypes, without loading the genotypes of each individual into R.
#'
#' Windows of `window_size` consecutive loci slide along the genome, one locus
#' at a time. A window is homozygous (for runs of homozygosity) if it contains
#' at most `max_opp_window` heterozygous genotypes and at most
#' `max_miss_window` missing genotypes, and no gap larger than `max_gap`
#' between consecutive loci. A locus is in a run if the proportion of
#' homozygous windows among all the windows that include it is larger than
#' `threshold`. Runs are stretches of consecutive loci in a run on a
#' chromosome, and only those with at least `min_snp` loci, a length of at
#' least `min_length_bps` and a density of at least `min_density` loci per
#' 1000 bps (and, if given, at most `max_opp_run` heterozygous and
#' `max_miss_run` missing genotypes) are returned. For runs of
#' heterozygosity, homozygous and heterozygous genotypes swap their roles.
#' Genotypes are classified by their dosage (0 and 2 are homozygous, 1 is
#' heterozygous).
#'
#' @param .x a [gen_tibble]
#' @param window_size the size of sliding window (number of SNP loci) (default =
#'   15)
#' @param threshold the threshold of overlapping windows of the same state
#'   (homozygous/heterozygous) to call a SNP in a RUN (default = 0.05)
#' @param min_snp minimum n. of SNP in a RUN (default = 3)
#' @param heterozygosity should we look for runs of heterozygosity (instead of
#'   homozygosity? (default = FALSE)
#' @param max_opp_window max n. of SNPs of the opposite type (e.g. heterozygous
#'   snps for runs of homozygosity) in the sliding window (default = 1)
#' @param max_miss_window max. n. of missing SNP in the sliding window (default
#'   = 1)
#' @param max_gap max distance between consecutive SNP to be still considered a
#'   potential run (default = 10^6 bps)
#' @param min_length_bps minimum length of run in bps (defaults to 1000 bps = 1
#'   kbps)
#' @param min_density minimum n. of SNP per kbps (defaults to 0.1 = 1 SNP every
#'   10 kbps)
#' @param max_opp_run max n. of opposite genotype SNPs in the run (optional)
#' @param max_miss_run max n. of missing SNPs in the run (optional)
#' @param n_cores the number of cores to use (individuals are analysed in
#'   parallel)
#'
#' @details This function returns a data frame with all runs detected in the
#'   dataset, in the format of `detectRUNS`, so that it can be used as the
#'   input of the functions of the `detectRUNS` package that create plots and
#'   produce statistics from the results (see plots and statistics functions
#'   in its manual, and/or refer to the `detectRUNS` vignette). The runs are
#'   the same as those of [detectRUNS::slidingRUNS.run()], including two
#'   features of its algorithm: windows slide across chromosomes (only the
#'   distance between consecutive loci in the loci table is checked against
#'   `max_gap`), and the counts of opposite and missing genotypes compared to
#'   `max_opp_run` and `max_miss_run` are those of the loci of the run. With
#'   fewer loci than `window_size`, no run is found. The loci must be sorted by
#'   chromosome and position (see [gt_order_loci()]); `detectRUNS` does not
#'   check this, and gives runs that start on the wrong chromosome when the
#'   loci of a chromosome are not adjacent.
#'
#'   If the [`gen_tibble`] is grouped, then the grouping variable is used to
#'   fill in the 'group' column (with several grouping variables, their values
#'   are pasted with ".", e.g. `pop1.F`). Otherwise, the 'group' column is
#'   filled with the same values as the 'id' column. Note that this behaviour is
#'   different from other windowed operations in `tidypopgen`, which
#'   return a list for grouped `gen_tibbles`; this different behaviour is
#'   designed to maintain compatibility with `detectRUNS`.
#'
#'   The old name for this function, `gt_roh_window`, is still available, but
#'   it is soft deprecated and will be removed in future versions of
#'   `tidypopgen`.
#'
#' @return A dataframe with RUNs of Homozygosity or Heterozygosity in the
#'   analysed dataset. The returned dataframe contains the following seven
#'   columns: "group", "id", "chrom", "nSNP", "from", "to", "lengthBps" (group:
#'   population, breed, case/control etc.; id: individual identifier; chrom:
#'   chromosome on which the run is located; nSNP: number of SNPs in the run;
#'   from: starting position of the run, in bps; to: end position of the run, in
#'   bps; lengthBps: size of the run). Runs are ordered by individual, then by
#'   chromosome (in the order of its levels) and position.
#' @export
#' @seealso [detectRUNS::slidingRUNS.run()], whose algorithm this function
#'   implements.
#' @references Biscarini F, Cozzi P, Gaspa G, Marras G (2019) `detectRUNS`:
#'   Detect Runs of Homozygosity and Runs of Heterozygosity in Diploid Genomes.
#'   R package. \doi{10.32614/CRAN.package.detectRUNS}
#' @examplesIf rlang::is_installed("detectRUNS")
#' sheep_ped <- system.file("extdata", "Kijas2016_Sheep_subset.ped",
#'   package = "detectRUNS"
#' )
#' sheep_gt <- tidypopgen::gen_tibble(sheep_ped,
#'   backingfile = tempfile(),
#'   quiet = TRUE
#' )
#' sheep_gt <- sheep_gt %>% group_by(population)
#' sheep_roh <- windows_indiv_roh(sheep_gt)
#' detectRUNS::plot_Runs(runs = sheep_roh)
windows_indiv_roh <- function(
  .x,
  window_size = 15,
  threshold = 0.05,
  min_snp = 3,
  heterozygosity = FALSE,
  max_opp_window = 1,
  max_miss_window = 1,
  max_gap = 10^6,
  min_length_bps = 1000,
  min_density = 1 / 1000,
  max_opp_run = NULL,
  max_miss_run = NULL,
  n_cores = bigstatsr::nb_cores()
) {
  stopifnot_has_genotypes(.x)
  if (!is.numeric(window_size) || length(window_size) != 1 ||
        is.na(window_size) || window_size < 1) {
    stop("window_size should be a positive integer")
  }
  loci <- show_loci(.x)
  chromosome <- loci$chromosome
  if (!is.factor(chromosome)) {
    chromosome <- as.factor(chromosome)
  }
  if (anyNA(chromosome)) {
    stop("the chromosomes of the loci should not be missing")
  }
  # windows slide along the loci table, so the loci of each chromosome must
  # be adjacent, and sorted by position (both checks are a single pass in
  # C++)
  positions_unsorted <- any_unsorted_by_chromosome(
    loci$position,
    chromosome
  )[["unsorted"]]
  if (chromosomes_not_adjacent(chromosome) || isTRUE(positions_unsorted)) {
    stop(
      "the loci should be sorted by chromosome and position ",
      "(see gt_order_loci())"
    )
  }
  # use groups (if defined)
  if (dplyr::is_grouped_df(.x)) {
    groups <- .group_of_rows(.x)
  } else {
    groups <- .x$id
  }

  # find the runs of all individuals in C++, in a single pass over the loci
  # (optional limits are NA when they are not used)
  runs <- roh_sliding_cpp(
    BM = .gt_get_fbm(.x),
    rowInd = .gt_fbm_rows(.x),
    colInd = .gt_fbm_cols(.x),
    chrom = chromosome,
    position = loci$position,
    window_size = window_size,
    threshold = threshold,
    heterozygosity = heterozygosity,
    max_opp_window = max_opp_window,
    max_miss_window = max_miss_window,
    max_gap = max_gap,
    min_snp = min_snp,
    min_length_bps = min_length_bps,
    min_density = min_density,
    max_opp_run = if (is.null(max_opp_run)) NA_real_ else max_opp_run,
    max_miss_run = if (is.null(max_miss_run)) NA_real_ else max_miss_run,
    ncores = n_cores,
    block_size = 10000L
  )

  # order the runs by individual, then by chromosome (in the order of its
  # levels) and position, as in detectRUNS
  run_order <- order(runs$indiv, runs$chrom, runs$start)
  indiv <- runs$indiv[run_order]
  from <- loci$position[runs$start[run_order]]
  to <- loci$position[runs$end[run_order]]
  data.frame(
    group = as.character(groups[indiv]),
    id = as.character(.x$id[indiv]),
    chrom = levels(chromosome)[runs$chrom[run_order]],
    nSNP = runs$n_snp[run_order],
    from = from,
    to = to,
    lengthBps = to - from,
    stringsAsFactors = FALSE
  )
}


# Alias for old name
#' @rdname windows_indiv_roh
#' @export
gt_roh_window <- function(
  .x,
  window_size = 15,
  threshold = 0.05,
  min_snp = 3,
  heterozygosity = FALSE,
  max_opp_window = 1,
  max_miss_window = 1,
  max_gap = 10^6,
  min_length_bps = 1000,
  min_density = 1 / 1000,
  max_opp_run = NULL,
  max_miss_run = NULL,
  n_cores = bigstatsr::nb_cores()
) {
  warning(
    "This is a soft-deprecated function, and will be removed in the ",
    "next version of tidypopgen. \n",
    "Please update your code to use windows_indiv_roh() instead. \n",
    "See ?windows_indiv_roh for more details. \n"
  )
  # call the new function
  windows_indiv_roh(
    .x = .x,
    window_size = window_size,
    threshold = threshold,
    min_snp = min_snp,
    heterozygosity = heterozygosity,
    max_opp_window = max_opp_window,
    max_miss_window = max_miss_window,
    max_gap = max_gap,
    min_length_bps = min_length_bps,
    min_density = min_density,
    max_opp_run = max_opp_run,
    max_miss_run = max_miss_run,
    n_cores = n_cores
  )
}
