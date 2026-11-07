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
#' @details The runs are the same as those of [detectRUNS::slidingRUNS.run()],
#'   including two features of its algorithm. First, windows slide along the
#'   loci table as a whole, so they cross from one chromosome to the next:
#'   near a boundary, a window holds the last loci of one chromosome and the
#'   first loci of the next. Only the distance between consecutive loci is
#'   checked against `max_gap`, and across a boundary it is usually negative
#'   (the next chromosome starts at a smaller position), so such windows are
#'   not discarded; they count towards the loci at the ends of both
#'   chromosomes, but runs never span two chromosomes, as each run ends with
#'   the last locus of its chromosome. Second, the counts of opposite and
#'   missing genotypes compared to `max_opp_run` and `max_miss_run` are those
#'   of the loci of the run. With fewer loci than `window_size`, no run is
#'   found. The loci must be sorted by chromosome and position (see
#'   [gt_order_loci()]).
#'
#'   Runs can be plotted with [autoplot.indiv_roh()], which draws the plots
#'   of `detectRUNS`. Use [as_detectRUNS()] to convert the runs to the format
#'   of `detectRUNS`, to use its own functions.
#'
#'   The old name for this function, `gt_roh_window`, is still available, but
#'   it is soft deprecated and will be removed in future versions of
#'   `tidypopgen`.
#'
#' @return An `indiv_roh` object, a tibble with one row per run, ordered by
#'   individual, then by chromosome (in the order of its levels) and
#'   position, with columns:
#'   - the grouping variable(s), if the [`gen_tibble`] is grouped (a single
#'   `group` column for one grouping variable, keeping its type, or one column
#'   per grouping variable);
#'   - `id`: the individual;
#'   - `chromosome`: the chromosome of the run (a factor, with the levels of
#'   the loci table);
#'   - `start`, `end`: the positions of the first and last locus of the run
#'   (in bp);
#'   - `length_bp`: `end - start`;
#'   - `n_loci`: the number of loci in the run;
#'   - `first_locus`, `last_locus`: the rows of the first and last locus of
#'   the run in the loci table ([show_loci()]);
#'   - `n_opp`: the number of loci of the run with the opposite genotype
#'   (heterozygous for runs of homozygosity, homozygous for runs of
#'   heterozygosity);
#'   - `n_missing`: the number of loci of the run with a missing genotype.
#'
#'   Its attributes `individuals` (a tibble of all the individuals, with
#'   their groups, including those with no run), `chromosome_length` (the
#'   length of each chromosome, its largest position, as in `detectRUNS`) and
#'   `heterozygosity` are used by [autoplot.indiv_roh()].
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
#' sheep_roh
#' # use the plots of detectRUNS
#' detectRUNS::plot_Runs(runs = as_detectRUNS(sheep_roh))
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
  first_locus <- runs$start[run_order]
  last_locus <- runs$end[run_order]
  start <- loci$position[first_locus]
  end <- loci$position[last_locus]
  res <- tibble::tibble(
    id = .x$id[indiv],
    chromosome = factor(levels(chromosome)[runs$chrom[run_order]],
      levels = levels(chromosome)
    ),
    start = start,
    end = end,
    length_bp = end - start,
    n_loci = runs$n_snp[run_order],
    first_locus = first_locus,
    last_locus = last_locus,
    n_opp = runs$n_opp[run_order],
    n_missing = runs$n_miss[run_order]
  )
  # the group of each run (as in the other tidy outputs of grouped
  # gen_tibbles), before the id; all individuals (with their groups) are
  # kept as an attribute, as some have no runs
  individuals <- tibble::tibble(id = .x$id)
  if (dplyr::is_grouped_df(.x)) {
    groups <- .group_info(.x)
    res <- dplyr::bind_cols(
      .group_columns(groups, groups$ids[indiv] + 1L),
      res
    )
    individuals <- dplyr::bind_cols(
      .group_columns(groups, groups$ids + 1L),
      individuals
    )
  }
  # the length of each chromosome, as in detectRUNS: its largest position,
  # i.e. the position of its last locus (loci are sorted); used for F_ROH
  chrom_runs <- vctrs::vec_unrep(chromosome)
  chromosome_length <- tibble::tibble(
    chromosome = chrom_runs$key,
    length_bp = loci$position[cumsum(chrom_runs$times)]
  )
  structure(res,
    individuals = individuals,
    chromosome_length = chromosome_length,
    heterozygosity = heterozygosity,
    class = c("indiv_roh", class(res))
  )
}


#' Convert runs of homozygosity to the format of `detectRUNS`
#'
#' Converts the runs found by [windows_indiv_roh()] to the data frame used by
#' the `detectRUNS` package, so that its functions to plot and summarise
#' runs (e.g. [detectRUNS::plot_Runs()], [detectRUNS::summaryRuns()]) can be
#' used.
#'
#' @param x a tibble of runs, as returned by [windows_indiv_roh()]
#' @returns a data frame with columns `group` (the grouping variables, pasted
#'   with "." if there are several, or the `id` if the runs were found on an
#'   ungrouped [`gen_tibble`]), `id`, `chrom`, `nSNP`, `from`, `to` and
#'   `lengthBps`, in the order of `x`
#' @export
#' @examplesIf rlang::is_installed("detectRUNS")
#' sheep_ped <- system.file("extdata", "Kijas2016_Sheep_subset.ped",
#'   package = "detectRUNS"
#' )
#' sheep_gt <- gen_tibble(sheep_ped, backingfile = tempfile(), quiet = TRUE)
#' sheep_roh <- windows_indiv_roh(sheep_gt %>% group_by(population))
#' head(as_detectRUNS(sheep_roh))
as_detectRUNS <- function(x) { # nolint: object_name_linter
  required <- c(
    "id", "chromosome", "start", "end", "length_bp", "n_loci"
  )
  if (!is.data.frame(x) || !all(required %in% names(x))) {
    stop(
      "x should be the output of windows_indiv_roh(), with columns ",
      paste(required, collapse = ", ")
    )
  }
  # the grouping variables are the columns before id
  group_cols <- names(x)[seq_len(match("id", names(x)) - 1)]
  if (length(group_cols) == 0) {
    group <- as.character(x$id)
  } else {
    group <- .group_labels(x[group_cols])
  }
  data.frame(
    group = group,
    id = as.character(x$id),
    chrom = as.character(x$chromosome),
    nSNP = x$n_loci,
    from = x$start,
    to = x$end,
    lengthBps = x$length_bp,
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
