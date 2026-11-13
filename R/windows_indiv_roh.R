#' Detect runs of homozygosity using a sliding-window
#'
#' This function looks for runs of homozygosity in a diploid genome with the
#' sliding-window algorithm of PLINK 1.9 (`--homozyg`), with the same
#' parameters and defaults. It is implemented natively in C++: all individuals
#' are analysed in a single pass over the genotypes.
#'
#' On each chromosome, a window of `window_snp` consecutive loci starts at
#' each locus (except the last `window_snp - 1`). A window is a hit if it has
#' at most `window_het` heterozygous and at most `window_missing` missing
#' genotypes. A locus is in a run if at least a proportion `window_threshold`
#' (rounded up) of the windows that include it are hits. Runs are stretches of
#' consecutive loci in a run, which are broken where the gap between two loci
#' is larger than `max_gap_kb`, and before the heterozygous genotype that
#' would exceed `max_het`. Runs are returned if they have at least `min_snp`
#' loci, a length of at least `min_kb`, and at most `max_kb_per_snp` kb per
#' locus.
#'
#' @param .x a [gen_tibble]
#' @param window_snp the number of loci in each window (PLINK
#'   `--homozyg-window-snp`, default 50)
#' @param window_het the largest number of heterozygous genotypes in a window
#'   that is a hit (PLINK `--homozyg-window-het`, default 1)
#' @param window_missing the largest number of missing genotypes in a window
#'   that is a hit (PLINK `--homozyg-window-missing`, default 5)
#' @param window_threshold the smallest proportion of the windows that
#'   include a locus that must be hits for the locus to be in a run (PLINK
#'   `--homozyg-window-threshold`, default 0.05)
#' @param min_snp the smallest number of loci in a run (PLINK `--homozyg-snp`,
#'   default 100)
#' @param min_kb the smallest length of a run, in kb (PLINK `--homozyg-kb`,
#'   default 1000)
#' @param max_kb_per_snp the largest average distance between the loci of a
#'   run, in kb per locus, i.e. the inverse of its smallest density (PLINK
#'   `--homozyg-density`, default 50)
#' @param max_gap_kb the largest gap between two consecutive loci of a run, in
#'   kb: a larger gap breaks a run (PLINK `--homozyg-gap`, default 1000)
#' @param max_het the largest number of heterozygous genotypes in a run
#'   (PLINK `--homozyg-het`, default `Inf`, i.e. no limit): a run is broken
#'   before the heterozygous genotype that would exceed it
#' @param subtract_1_from_lengths if `TRUE`, the length of a run is the
#'   distance between its first and last positions; by default it is that
#'   distance + 1, as in PLINK 1.9 (PLINK `--homozyg subtract-1-from-lengths`,
#'   which gives the lengths of PLINK 1.07)
#' @param n_cores the number of cores to use (individuals are analysed in
#'   parallel)
#' @param ... not used. The arguments of the version of this function based
#'   on `detectRUNS` (e.g. `window_size`) give an error: use
#'   [gt_roh_window()] for that version.
#'
#' @details The runs are those of `plink --homozyg` with the same parameters,
#'   including the features of its algorithm:
#'   - windows restart at each chromosome, so a chromosome with fewer loci
#'   than `window_snp` has no runs; windows can span large gaps, which only
#'   break runs;
#'   - the number of windows used for the threshold of the last
#'   `window_snp - 1` loci of a chromosome is the number of loci to its end
#'   (which, on chromosomes with fewer than `2 * window_snp - 2` loci, is
#'   more than the number of windows that include the locus);
#'   - runs are not trimmed, so they can start or end with heterozygous or
#'   missing genotypes; missing genotypes in a run are not limited.
#'
#'   Genotypes are classified by their dosage (1 is heterozygous, 0 and 2 are
#'   homozygous). Unlike PLINK, all chromosomes are analysed: PLINK skips
#'   haploid chromosomes (e.g. Y and MT), and male X chromosomes. Remove them
#'   first (e.g. with [select_loci()]) to get the same runs. PLINK has no
#'   runs of heterozygosity; [gt_roh_window()] (based on `detectRUNS`) can
#'   still find them. The loci must be sorted by chromosome and position (see
#'   [gt_order_loci()]).
#'
#'   Runs can be plotted with [autoplot.indiv_roh()], which draws the plots
#'   of `detectRUNS`. Use [as_detectRUNS()] to convert the runs to the format
#'   of `detectRUNS`, to use its own functions.
#'
#' @return An `indiv_roh` object, a tibble with one row per run, ordered by
#'   individual, then by chromosome (in the order of its levels) and
#'   position, with columns:
#'   - the grouping variable(s), if the [`gen_tibble`] is grouped (a single
#'   `group` column for one grouping variable, keeping its type, or one column
#'   per grouping variable);
#'   - `id`: the individual;
#'   - `chromosome`: the chromosome of the run (a factor, with the levels of
#'   the loci table; integer chromosomes, from `chromosomes_as_int = TRUE`,
#'   are replaced by their names from the chromosome sidecar, see
#'   [gt_lookup_chromosomes()], so the runs are the same as without
#'   `chromosomes_as_int`);
#'   - `start`, `end`: the positions of the first and last locus of the run
#'   (in bp; PLINK `POS1` and `POS2`);
#'   - `length_bp`: `end - start + 1` (or `end - start` with
#'   `subtract_1_from_lengths = TRUE`; PLINK `KB` x 1000);
#'   - `n_loci`: the number of loci in the run (PLINK `NSNP`);
#'   - `first_locus`, `last_locus`: the rows of the first and last locus of
#'   the run in the loci table ([show_loci()]);
#'   - `n_het`: the number of heterozygous genotypes in the run (PLINK `PHET`
#'   x `NSNP`);
#'   - `n_missing`: the number of missing genotypes in the run.
#'
#'   Its attributes `individuals` (a tibble of all the individuals, with
#'   their groups, including those with no run), `chromosome_length` (the
#'   length of each chromosome, its largest position, as in `detectRUNS`) and
#'   `n_loci` (the number of loci) are used by [autoplot.indiv_roh()],
#'   [gt_roh_froh()], [gt_roh_loci()] and [gt_roh_islands()].
#' @export
#' @references Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee
#'   JJ (2015) Second-generation PLINK: rising to the challenge of larger and
#'   richer datasets. GigaScience, 4. \doi{10.1186/s13742-015-0047-8}
#' @examples
#' # a simulated dataset of three populations with increasing inbreeding
#' # (see data-raw/create_roh_dataset.R in the source of the package)
#' bed_file <- system.file("extdata", "roh", "roh_sim.bed",
#'   package = "tidypopgen"
#' )
#' sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
#' sim_gt <- sim_gt %>% group_by(population)
#' sim_roh <- windows_indiv_roh(sim_gt)
#' sim_roh
#' # the number of runs in each population
#' table(sim_roh$group)
#' \donttest{
#' # use the plots of detectRUNS
#' if (rlang::is_installed("detectRUNS")) {
#'   detectRUNS::plot_Runs(runs = as_detectRUNS(sim_roh))
#' }
#' }
windows_indiv_roh <- function(
  .x,
  window_snp = 50,
  window_het = 1,
  window_missing = 5,
  window_threshold = 0.05,
  min_snp = 100,
  min_kb = 1000,
  max_kb_per_snp = 50,
  max_gap_kb = 1000,
  max_het = Inf,
  subtract_1_from_lengths = FALSE,
  n_cores = bigstatsr::nb_cores(),
  ...
) {
  # the arguments of the old version (based on detectRUNS) are not used
  # anymore, as the parameters are now those of PLINK
  old_args <- intersect(
    names(list(...)),
    c(
      "window_size", "threshold", "heterozygosity", "max_opp_window",
      "max_miss_window", "max_gap", "min_length_bps", "min_density",
      "max_opp_run", "max_miss_run"
    )
  )
  if (length(old_args) > 0) {
    stop(
      "windows_indiv_roh() now follows PLINK (--homozyg), and argument(s) ",
      paste(old_args, collapse = ", "), " are no longer used: see ",
      "?windows_indiv_roh for the new ones, or use gt_roh_window() for the ",
      "old version, based on detectRUNS",
      call. = FALSE
    )
  }
  rlang::check_dots_empty()
  stopifnot_has_genotypes(.x)
  is_count <- function(x, min = 0) {
    is.numeric(x) && length(x) == 1 && !is.na(x) && x >= min
  }
  if (!is_count(window_snp, 1)) {
    stop("window_snp should be a positive integer")
  }
  for (arg in c("window_het", "window_missing", "min_snp", "min_kb",
                "max_kb_per_snp", "max_gap_kb", "max_het")) {
    if (!is_count(get(arg))) {
      stop(arg, " should be a non-negative number")
    }
  }
  if (!is_count(window_threshold) || window_threshold > 1) {
    stop("window_threshold should be between 0 and 1")
  }
  loci <- show_loci(.x)
  # integer chromosomes (chromosomes_as_int) with a chromosome sidecar become
  # the factor that a read without chromosomes_as_int gives, so that the runs
  # are ordered by chromosome name and the outputs hold names (only the
  # distinct codes are looked up, and the factor needs no string per locus)
  chromosome <- roh_chromosome_factor(.x, loci$chromosome)
  # integer chromosomes without a sidecar (e.g. numbers in a loci table) are
  # used as codes (a factor would create a string for each chromosome)
  int_chromosomes <- is.integer(chromosome) && !is.factor(chromosome)
  if (!int_chromosomes && !is.factor(chromosome)) {
    chromosome <- as.factor(chromosome)
  }
  # (the C++ code takes positive codes, with a table as long as the largest:
  # integer chromosomes are turned into groups 1..n if their codes are large)
  chromosome_groups <- chromosome_group_codes(chromosome)
  if (anyNA(chromosome)) {
    stop("the chromosomes of the loci should not be missing")
  }
  if (anyNA(loci$position)) {
    stop(
      sum(is.na(loci$position)), " loci have a missing position, ",
      "so runs can not be found"
    )
  }
  # windows slide along each chromosome in a single pass over the loci
  # table, so the loci of each chromosome must be adjacent, and sorted by
  # position (both checks are a single pass in C++)
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

  # find the runs of all individuals in C++, in a single pass over the loci;
  # kb are converted to bp as in PLINK (integers for the length and the gap,
  # with a small tolerance for the density)
  plink_epsilon <- 2^-30
  runs <- roh_sliding_cpp(
    BM = .gt_get_fbm(.x),
    rowInd = .gt_fbm_rows(.x),
    colInd = .gt_fbm_cols(.x),
    chrom = chromosome_groups$codes,
    position = loci$position,
    window_size = as.integer(window_snp),
    max_het_window = as.integer(min(window_het, .Machine$integer.max)),
    max_miss_window = as.integer(min(window_missing, .Machine$integer.max)),
    threshold = window_threshold,
    max_gap = floor(max_gap_kb * 1000 + plink_epsilon),
    min_snp = min_snp,
    min_bases = floor(min_kb * 1000 + plink_epsilon),
    max_bases_per_snp = max_kb_per_snp * 1000 + plink_epsilon,
    max_het = max_het,
    old_lengths = isTRUE(subtract_1_from_lengths),
    ncores = n_cores,
    block_size = 10000L
  )

  # order the runs by individual, then by chromosome (in the order of its
  # levels) and position
  run_order <- order(runs$indiv, runs$chrom, runs$start)
  indiv <- runs$indiv[run_order]
  first_locus <- runs$start[run_order]
  last_locus <- runs$end[run_order]
  start <- loci$position[first_locus]
  end <- loci$position[last_locus]
  res <- tibble::tibble(
    id = .x$id[indiv],
    chromosome = if (int_chromosomes) {
      if (is.null(chromosome_groups$values)) {
        runs$chrom[run_order]
      } else {
        chromosome_groups$values[runs$chrom[run_order]]
      }
    } else {
      factor(levels(chromosome)[runs$chrom[run_order]],
        levels = levels(chromosome)
      )
    },
    start = start,
    end = end,
    length_bp = end - start + !isTRUE(subtract_1_from_lengths),
    n_loci = runs$n_snp[run_order],
    first_locus = first_locus,
    last_locus = last_locus,
    n_het = runs$n_het[run_order],
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
    n_loci = nrow(loci),
    class = c("indiv_roh", class(res))
  )
}


# the chromosomes of the loci of .x for the runs of homozygosity: integer
# chromosomes (chromosomes_as_int) become a factor of their names, from the
# chromosome sidecar, with the levels that a read without chromosomes_as_int
# gives (see chromosome_sidecar_factor()), so that the runs are the same
# with and without integer chromosomes. Other chromosomes, and integer
# chromosomes without a sidecar (e.g. numbers from a bed file), are returned
# unchanged
roh_chromosome_factor <- function(.x, chromosome) {
  if (is_int_column(chromosome) &&
        file.exists(chromosomes_sidecar_path(.x))) {
    chromosome <- chromosome_sidecar_factor(.x, chromosome)
  }
  chromosome
}

# the loci table given to gt_roh_loci() and gt_roh_islands(), with integer
# chromosomes (chromosomes_as_int) replaced by the chromosomes of the runs,
# if windows_indiv_roh() found them as names (see roh_chromosome_factor()):
# the attribute chromosome_length of roh has a row for each stretch of
# loci of the same chromosome, in the order of the loci table, so each
# stretch of loci gets the chromosome of its row (after checking that the
# stretches and their last positions match)
roh_loci_chromosomes <- function(roh, loci) {
  chromosome_length <- attr(roh, "chromosome_length")
  if (!is_int_column(loci$chromosome) ||
        is_int_column(chromosome_length$chromosome)) {
    return(loci)
  }
  chrom_runs <- vctrs::vec_unrep(loci$chromosome)
  last_positions <- as.numeric(loci$position[cumsum(chrom_runs$times)])
  if (nrow(chrom_runs) != nrow(chromosome_length) ||
        !identical(last_positions, as.numeric(chromosome_length$length_bp))) {
    stop(
      "loci should be the loci table (show_loci()) of the gen_tibble ",
      "used to find the runs"
    )
  }
  loci$chromosome <- rep(chromosome_length$chromosome, times = chrom_runs$times)
  loci
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
#' @examples
#' bed_file <- system.file("extdata", "roh", "roh_sim.bed",
#'   package = "tidypopgen"
#' )
#' sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
#' sim_roh <- windows_indiv_roh(sim_gt %>% group_by(population))
#' head(as_detectRUNS(sim_roh))
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


#' Detect runs of homozygosity with `detectRUNS` (deprecated)
#'
#' This function is deprecated: use [windows_indiv_roh()], which finds runs
#' of homozygosity as PLINK (`--homozyg`), with its parameters. It is kept so
#' that old code still gives the same results: it calls the sliding-window
#' algorithm of the package `detectRUNS` ([detectRUNS::slidingRUNS.run()]) on
#' each individual, as `windows_indiv_roh()` used to do. It needs `detectRUNS`
#' to be installed, and will be removed in a future version of `tidypopgen`.
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
#'
#' @details If the [`gen_tibble`] is grouped, then the grouping variable is
#'   used to fill in the 'group' column (with several grouping variables,
#'   their values are pasted with ".", e.g. `pop1.F`). Otherwise, the 'group'
#'   column is filled with the same values as the 'id' column.
#'
#' @return A dataframe with RUNs of Homozygosity or Heterozygosity in the
#'   analysed dataset, in the format of `detectRUNS`, with columns "group",
#'   "id", "chrom", "nSNP", "from", "to", "lengthBps" (group: population,
#'   breed, case/control etc.; id: individual identifier; chrom: chromosome on
#'   which the run is located; nSNP: number of SNPs in the run; from: starting
#'   position of the run, in bps; to: end position of the run, in bps;
#'   lengthBps: size of the run)
#' @export
#' @keywords internal
#' @seealso [detectRUNS::slidingRUNS.run()] which this function wraps.
#' @references Biscarini F, Cozzi P, Gaspa G, Marras G (2019) `detectRUNS`:
#'   Detect Runs of Homozygosity and Runs of Heterozygosity in Diploid Genomes.
#'   R package. \doi{10.32614/CRAN.package.detectRUNS}
#' @examplesIf rlang::is_installed("detectRUNS")
#' bed_file <- system.file("extdata", "roh", "roh_sim.bed",
#'   package = "tidypopgen"
#' )
#' sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
#' # deprecated: use windows_indiv_roh() instead
#' sim_roh <- suppressWarnings(gt_roh_window(sim_gt %>% group_by(population)))
#' head(sim_roh)
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
  max_miss_run = NULL
) {
  warning(
    "gt_roh_window() is deprecated, and will be removed in a future ",
    "version of tidypopgen. \n",
    "Please use windows_indiv_roh(), which follows PLINK (--homozyg), ",
    "instead. \n",
    "See ?windows_indiv_roh for more details. \n"
  )
  if (!requireNamespace("detectRUNS", quietly = TRUE)) {
    stop(
      "to use this function, first install package 'detectRUNS' with\n",
      "install.packages('detectRUNS')"
    )
  }
  # collect all parameters in a variable
  parameters <- list(
    windowSize = window_size,
    threshold = threshold,
    minSNP = min_snp,
    ROHet = heterozygosity,
    maxOppWindow = max_opp_window,
    maxMissWindow = max_miss_window,
    maxGap = max_gap,
    minLengthBps = min_length_bps,
    minDensity = min_density,
    maxOppRun = max_opp_run,
    maxMissRun = max_miss_run
  )

  # create a map object
  map <- show_loci(.x) %>%
    dplyr::select(dplyr::all_of(c("chromosome", "name", "position"))) %>%
    dplyr::rename(
      "Chrom" = "chromosome",
      "SNP" = "name",
      "bps" = "position"
    ) %>%
    as.data.frame()
  # integer chromosomes (chromosomes_as_int) get their names, as without
  # chromosomes_as_int
  map$Chrom <- roh_chromosome_factor(.x, map$Chrom)

  # compute the gaps between snps
  gaps <- diff(map$bps)
  # use groups (if defined)
  if (dplyr::is_grouped_df(.x)) {
    groups <- .group_of_rows(.x)
  } else {
    groups <- .x$id
  }

  runs_list <- vector("list", nrow(.x))
  X <- .gt_get_fbm(.x) # pointer for the FBM #nolint
  col_ind <- .gt_fbm_cols(.x) # column indices for the snps to consider
  rows_to_keep <- .gt_fbm_rows(.x)
  for (i in seq_len(nrow(.x))) {
    this_genotype <- X[rows_to_keep[i], col_ind]
    # detectRUNS expects 0 for homozygotes and 1 for heterozygotes (its
    # readers recode 0/1/2 dosages with genoConvert), otherwise homozygotes
    # for the alternate allele count as heterozygotes for runs of
    # heterozygosity
    this_genotype[which(this_genotype == 2)] <- 0
    this_indiv <- list(FID = groups[i], IID = .x$id[i])
    # find runs for this individual
    runs_list[[i]] <-
      utils::getFromNamespace("slidingRuns", "detectRUNS")(
        this_genotype,
        this_indiv,
        map,
        gaps,
        parameters
      ) # nolint
  }

  # return calculated runs (data.frame)
  runs_df <- dplyr::bind_rows(runs_list)
  # fix row names
  row.names(runs_df) <- NULL
  runs_df
}
