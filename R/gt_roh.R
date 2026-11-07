#' Inbreeding from runs of homozygosity (F_ROH)
#'
#' Computes F_ROH, the proportion of the genome of each individual that is in
#' runs of homozygosity, from the runs found by [windows_indiv_roh()]. F_ROH
#' can be computed from the runs of at least a given length (e.g. 1, 2, 4, 8
#' and 16 Mb): long runs come from recent common ancestors, while short runs
#' also reflect older inbreeding.
#'
#' The length of the genome is, by default, the sum of the lengths of the
#' chromosomes, each being its largest position (as in `detectRUNS`); it can
#' be given instead with `genome_length` (e.g. the length of the autosomes of
#' the reference genome, or the length covered by the loci).
#'
#' @param roh an `indiv_roh` object, as returned by [windows_indiv_roh()]
#' @param min_length_mb the smallest length of the runs to use, in Mb (default
#'   0, i.e. all runs). With several values, F_ROH is computed for each.
#' @param by whether to compute F_ROH over the whole genome (`genome`, the
#'   default), or for each chromosome (`chromosome`)
#' @param genome_length the length of the genome, in bp (only used with `by
#'   = "genome"`). Default is the sum of the lengths of the chromosomes (see
#'   Details).
#' @returns a tibble with one row per individual (including those with no
#'   runs), with the grouping variables of `roh` (if any), `id`, the column
#'   `min_length_mb` if there are several values, `chromosome` if `by =
#'   "chromosome"`, `length_bp` (the total length of the runs) and `froh`.
#' @export
#' @seealso [autoplot.indiv_roh()] to plot F_ROH.
#' @examples
#' bed_file <- system.file("extdata", "roh", "roh_sim.bed",
#'   package = "tidypopgen"
#' )
#' sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
#' sim_roh <- windows_indiv_roh(sim_gt %>% group_by(population))
#' # F_ROH of each individual
#' gt_roh_froh(sim_roh)
#' # from the runs of at least 1, 2 and 4 Mb
#' gt_roh_froh(sim_roh, min_length_mb = c(1, 2, 4))
#' # on each chromosome
#' gt_roh_froh(sim_roh, by = "chromosome")
gt_roh_froh <- function(roh, min_length_mb = 0,
                        by = c("genome", "chromosome"),
                        genome_length = NULL) {
  by <- match.arg(by)
  .check_indiv_roh(roh)
  if (!is.numeric(min_length_mb) || length(min_length_mb) < 1 ||
        anyNA(min_length_mb)) {
    stop("min_length_mb should be a numeric vector")
  }
  chromosome_length <- attr(roh, "chromosome_length")
  if (by == "genome") {
    if (is.null(genome_length)) {
      genome_length <- sum(as.numeric(chromosome_length$length_bp))
    } else if (!is.numeric(genome_length) || length(genome_length) != 1 ||
                 is.na(genome_length) || genome_length <= 0) {
      stop("genome_length should be a positive number")
    }
  }
  individuals <- attr(roh, "individuals")
  runs <- tibble::as_tibble(roh)
  res <- lapply(min_length_mb, function(min_length) {
    long_runs <- runs[runs$length_bp >= min_length * 1e6, ]
    if (by == "genome") {
      # the total length of the runs of each individual (0 if none)
      total <- tapply(
        as.numeric(long_runs$length_bp),
        factor(long_runs$id, levels = individuals$id),
        sum,
        default = 0
      )
      out <- individuals
      out$length_bp <- as.numeric(total)
      out$froh <- out$length_bp / genome_length
    } else {
      # the total length of the runs of each individual on each chromosome
      chromosomes <- chromosome_length$chromosome
      total <- tapply(
        as.numeric(long_runs$length_bp),
        list(
          factor(long_runs$id, levels = individuals$id),
          factor(as.character(long_runs$chromosome),
            levels = as.character(chromosomes)
          )
        ),
        sum,
        default = 0
      )
      rows <- rep(seq_len(nrow(individuals)), times = length(chromosomes))
      out <- individuals[rows, ]
      out$chromosome <- rep(chromosomes, each = nrow(individuals))
      out$length_bp <- as.numeric(total)
      out$froh <- out$length_bp /
        rep(as.numeric(chromosome_length$length_bp), each = nrow(individuals))
      out <- out[order(match(out$id, individuals$id)), ]
    }
    if (length(min_length_mb) > 1) {
      out <- tibble::add_column(out, min_length_mb = min_length,
        .after = "id"
      )
    }
    out
  })
  res <- dplyr::bind_rows(res)
  if (length(min_length_mb) > 1) {
    res <- res[order(match(res$id, individuals$id)), ]
  }
  res
}

#' The loci in runs of homozygosity
#'
#' Counts, for each locus, the individuals that are in a run of homozygosity
#' (as found by [windows_indiv_roh()]), overall or in each group. Loci that
#' are in runs in most individuals of a population (see [gt_roh_islands()])
#' can point to regions under selection. The counts are those of the
#' `.hom.summary` file of PLINK `--homozyg`.
#'
#' @param roh an `indiv_roh` object, as returned by [windows_indiv_roh()]
#' @param loci the loci table of the [`gen_tibble`] used to find the runs
#'   (i.e. `show_loci(x)`, as the runs refer to its rows)
#' @param by_group if `roh` was found on a grouped [`gen_tibble`], whether to
#'   count the individuals of each group (the default) or of all groups
#' @returns a tibble with one row per locus (and per group), with the grouping
#'   variables (if counting by group), `name`, `chromosome`, `position`,
#'   `n_roh` (the number of individuals in a run), `n_indiv` (the number of
#'   individuals) and `prop_roh` (`n_roh / n_indiv`)
#' @export
#' @seealso [gt_roh_islands()] to find the regions in runs in many
#'   individuals, and [autoplot.indiv_roh()] (`type = "coverage"`) to plot
#'   them.
#' @examples
#' bed_file <- system.file("extdata", "roh", "roh_sim.bed",
#'   package = "tidypopgen"
#' )
#' sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
#' sim_roh <- windows_indiv_roh(sim_gt %>% group_by(population))
#' roh_loci <- gt_roh_loci(sim_roh, show_loci(sim_gt))
#' # the loci most often in runs in each population
#' roh_loci %>%
#'   group_by(group) %>%
#'   slice_max(prop_roh, n = 1, with_ties = FALSE)
gt_roh_loci <- function(roh, loci, by_group = TRUE) {
  .check_indiv_roh(roh)
  .check_roh_loci(roh, loci)
  n_loci <- nrow(loci)
  individuals <- attr(roh, "individuals")
  group_cols <- setdiff(names(individuals), "id")
  if (!by_group) group_cols <- character(0)
  # the groups (as rows of a table of the grouping variables), and the group
  # of each run and of each individual
  if (length(group_cols) > 0) {
    keys <- dplyr::distinct(individuals[group_cols])
    keys <- dplyr::arrange(keys, dplyr::across(dplyr::everything()))
    indiv_group <- vctrs::vec_match(individuals[group_cols], keys)
    run_group <- indiv_group[match(roh$id, individuals$id)]
  } else {
    keys <- NULL
    indiv_group <- rep(1L, nrow(individuals))
    run_group <- rep(1L, nrow(roh))
  }
  n_groups <- max(c(indiv_group, 1L))
  res <- lapply(seq_len(n_groups), function(g) {
    in_group <- run_group == g
    # each run adds 1 from its first locus to its last (a difference
    # array, summed along the loci)
    n_roh <- cumsum(
      tabulate(roh$first_locus[in_group], nbins = n_loci + 1) -
        tabulate(roh$last_locus[in_group] + 1L, nbins = n_loci + 1)
    )[seq_len(n_loci)]
    n_indiv <- sum(indiv_group == g)
    out <- tibble::tibble(
      name = loci$name,
      chromosome = loci$chromosome,
      position = loci$position,
      n_roh = as.integer(n_roh),
      n_indiv = n_indiv,
      prop_roh = n_roh / n_indiv
    )
    if (!is.null(keys)) {
      out <- dplyr::bind_cols(keys[rep(g, n_loci), ], out)
    }
    out
  })
  dplyr::bind_rows(res)
}

#' Islands of runs of homozygosity
#'
#' Finds the regions of the genome where many individuals (at least a
#' proportion `threshold` of them, overall or in each group) are in runs of
#' homozygosity: islands of runs of homozygosity (also called ROH hotspots),
#' which can point to regions under selection.
#'
#' Islands are stretches of consecutive loci of a chromosome that are in runs
#' in at least `threshold` of the individuals (see [gt_roh_loci()]).
#'
#' @inheritParams gt_roh_loci
#' @param threshold the smallest proportion of individuals in a run for a
#'   locus to be in an island (default 0.5)
#' @returns a tibble with one row per island, with the grouping variables (if
#'   counting by group), `chromosome`, `start` and `end` (the positions of its
#'   first and last loci), `length_bp` (`end - start + 1`), `n_loci`,
#'   `first_locus` and `last_locus` (rows of the loci table), and `max_prop`
#'   and `mean_prop` (the largest and mean proportions of individuals in a
#'   run over its loci)
#' @export
#' @examples
#' bed_file <- system.file("extdata", "roh", "roh_sim.bed",
#'   package = "tidypopgen"
#' )
#' sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
#' sim_roh <- windows_indiv_roh(sim_gt %>% group_by(population))
#' # regions in runs in at least 20% of the individuals of a population
#' gt_roh_islands(sim_roh, show_loci(sim_gt), threshold = 0.2)
gt_roh_islands <- function(roh, loci, threshold = 0.5, by_group = TRUE) {
  if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) ||
        threshold <= 0 || threshold > 1) {
    stop("threshold should be a proportion larger than 0")
  }
  roh_loci <- gt_roh_loci(roh, loci, by_group = by_group)
  n_loci <- nrow(loci)
  group_cols <- names(roh_loci)[seq_len(match("name", names(roh_loci)) - 1)]
  n_groups <- nrow(roh_loci) / n_loci
  res <- lapply(seq_len(n_groups), function(g) {
    rows <- (g - 1) * n_loci + seq_len(n_loci)
    prop <- roh_loci$prop_roh[rows]
    in_island <- which(prop >= threshold)
    if (length(in_island) == 0) {
      return(NULL)
    }
    # a new island starts after a locus that is not in an island, or on a
    # new chromosome
    chrom <- as.integer(loci$chromosome)
    gap <- diff(in_island) != 1 | diff(chrom[in_island]) != 0
    new_island <- c(TRUE, gap)
    island <- cumsum(new_island)
    first <- in_island[new_island]
    last <- in_island[c(new_island[-1], TRUE)]
    out <- tibble::tibble(
      chromosome = loci$chromosome[first],
      start = loci$position[first],
      end = loci$position[last],
      length_bp = loci$position[last] - loci$position[first] + 1L,
      n_loci = as.integer(last - first + 1L),
      first_locus = first,
      last_locus = last,
      max_prop = as.numeric(tapply(prop[in_island], island, max)),
      mean_prop = as.numeric(tapply(prop[in_island], island, mean))
    )
    if (length(group_cols) > 0) {
      out <- dplyr::bind_cols(
        roh_loci[rep(rows[1], nrow(out)), group_cols],
        out
      )
    }
    out
  })
  res <- dplyr::bind_rows(res)
  if (nrow(res) == 0) {
    # no island: an empty tibble with the columns of the islands
    res <- tibble::tibble(
      chromosome = loci$chromosome[0], start = loci$position[0],
      end = loci$position[0], length_bp = loci$position[0],
      n_loci = integer(0), first_locus = integer(0),
      last_locus = integer(0), max_prop = numeric(0),
      mean_prop = numeric(0)
    )
    if (length(group_cols) > 0) {
      res <- dplyr::bind_cols(roh_loci[0, group_cols], res)
    }
  }
  res
}

# check that roh is the output of windows_indiv_roh(), with its attributes
.check_indiv_roh <- function(roh) {
  if (!inherits(roh, "indiv_roh") ||
        is.null(attr(roh, "individuals")) ||
        is.null(attr(roh, "chromosome_length"))) {
    stop(
      "roh should be the output of windows_indiv_roh() (some operations ",
      "on it drop the information on the individuals and chromosomes)"
    )
  }
}

# check that loci is the loci table of the gen_tibble of the runs
.check_roh_loci <- function(roh, loci) {
  if (!is.data.frame(loci) ||
        !all(c("name", "chromosome", "position") %in% names(loci))) {
    stop("loci should be the loci table (show_loci()) of the gen_tibble")
  }
  n_loci <- attr(roh, "n_loci")
  if ((!is.null(n_loci) && nrow(loci) != n_loci) ||
        (nrow(roh) > 0 && max(roh$last_locus) > nrow(loci))) {
    stop(
      "loci should be the loci table (show_loci()) of the gen_tibble ",
      "used to find the runs"
    )
  }
}
