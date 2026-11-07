#' Autoplots for `indiv_roh` objects
#'
#' Plots of runs of homozygosity (or heterozygosity) found by
#' [windows_indiv_roh()], equivalent to the plots of `detectRUNS`. Runs are
#' coloured (or split) by group if the [`gen_tibble`] was grouped. The
#' following types of plots are available:
#' - `runs`: the runs of each individual along each chromosome (as
#'   [detectRUNS::plot_Runs()]);
#' - `stacked`: the runs of each group stacked along each chromosome, one run
#'   per line, ordered by start (as [detectRUNS::plot_StackedRuns()]);
#' - `coverage`: the percentage of individuals of each group in a run along
#'   each chromosome (as [detectRUNS::plot_SnpsInRuns()] and
#'   [detectRUNS::plot_manhattanRuns()], computed from the runs rather than
#'   at each locus);
#' - `froh`: the inbreeding coefficient F_ROH, the total length of the runs
#'   of an individual over the length of the genome (with `froh_by =
#'   "individual"`, one point per individual, by group) or of each chromosome
#'   (with `froh_by = "chromosome"`, the mean of each group) (as
#'   [detectRUNS::plot_InbreedingChr()]); the length of a chromosome is its
#'   largest position, as in `detectRUNS`;
#' - `violin`: the distribution of the total (or mean, see `method`) length of
#'   the runs of each individual, by group (as [detectRUNS::plot_ViolinRuns()]);
#' - `length_class`: the percentage of runs in each class of length, by group
#'   (as [detectRUNS::plot_DistributionRuns()]).
#'
#' Individuals with no runs count in `coverage`, `froh` and `violin` (with no
#' runs, or a length of 0).
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are not
#' customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' @param object an object of class `indiv_roh`, as returned by
#'   [windows_indiv_roh()]
#' @param type the type of plot (`runs`, `stacked`, `coverage`, `froh`,
#'   `violin` or `length_class`)
#' @param chromosome the chromosomes to plot (by default, all those with runs;
#'   for `froh` with `froh_by = "individual"`, the genome is always used)
#' @param froh_by for `froh`, whether to plot F_ROH of each individual over
#'   the genome (`individual`), or the mean F_ROH of each group on each
#'   chromosome (`chromosome`)
#' @param method for `violin`, whether to plot the total (`sum`) or the mean
#'   (`mean`) length of the runs of each individual
#' @param classes for `length_class`, the lower bounds of the classes of length
#'   of the runs, in Mb (the last class has no upper bound)
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @export
#' @examples
#' # a simulated dataset of three populations with increasing inbreeding
#' bed_file <- system.file("extdata", "roh", "roh_sim.bed",
#'   package = "tidypopgen"
#' )
#' sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
#' sim_roh <- windows_indiv_roh(sim_gt %>% group_by(population))
#'
#' autoplot(sim_roh, type = "runs", chromosome = "3")
#' autoplot(sim_roh, type = "coverage")
#' autoplot(sim_roh, type = "froh")
#' autoplot(sim_roh, type = "length_class")
autoplot.indiv_roh <- function(
  object,
  type = c("runs", "stacked", "coverage", "froh", "violin", "length_class"),
  chromosome = NULL,
  froh_by = c("individual", "chromosome"),
  method = c("sum", "mean"),
  classes = c(0, 2, 4, 8, 16),
  ...
) {
  rlang::check_dots_empty()
  type <- match.arg(type)
  froh_by <- match.arg(froh_by)
  method <- match.arg(method)

  runs <- .roh_runs_with_group(object)
  individuals <- .roh_individuals(object)
  chromosome_length <- attr(object, "chromosome_length")
  if (is.null(chromosome_length)) {
    # (attributes are lost by some operations on the runs; the largest end
    # of a run is then used as the length of each chromosome)
    chromosome_length <- runs %>%
      dplyr::group_by(.data$chromosome) %>%
      dplyr::summarise(length_bp = max(.data$end))
  }
  # the chromosomes to plot
  if (!is.null(chromosome)) {
    if (!all(chromosome %in% as.character(runs$chromosome))) {
      stop("chromosome should be among the chromosomes with runs")
    }
    runs <- runs[as.character(runs$chromosome) %in% chromosome, ]
    if (!(type == "froh" && froh_by == "individual")) {
      chromosome_length <- chromosome_length[
        as.character(chromosome_length$chromosome) %in% chromosome,
      ]
    }
  }
  run_word <- if (isTRUE(attr(object, "heterozygosity"))) {
    "Runs of heterozygosity"
  } else {
    "Runs of homozygosity"
  }

  if (type == "runs") {
    # one line per individual (ordered by group), one segment per run
    runs$id <- factor(as.character(runs$id), levels = individuals$id)
    p <- ggplot2::ggplot(runs) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = .data$start / 1e6, xend = .data$end / 1e6,
          y = .data$id, yend = .data$id, colour = .data$group
        ),
        linewidth = 2
      ) +
      ggplot2::facet_wrap(~ chromosome, scales = "free_x") +
      ggplot2::labs(
        x = "Position (Mb)", y = "Individual", colour = "Group",
        title = run_word
      )
  } else if (type == "stacked") {
    # each run on its own line, ordered by start, within each chromosome
    # and group
    runs <- runs %>%
      dplyr::arrange(.data$chromosome, .data$group, .data$start) %>%
      dplyr::group_by(.data$chromosome, .data$group) %>%
      dplyr::mutate(line = dplyr::row_number()) %>%
      dplyr::ungroup()
    p <- ggplot2::ggplot(runs) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = .data$start / 1e6, xend = .data$end / 1e6,
          y = .data$line, yend = .data$line, colour = .data$group
        )
      ) +
      ggplot2::facet_grid(group ~ chromosome, scales = "free") +
      ggplot2::labs(
        x = "Position (Mb)", y = "Runs", colour = "Group",
        title = run_word
      )
  } else if (type == "coverage") {
    coverage <- .roh_coverage(runs, individuals)
    p <- ggplot2::ggplot(coverage) +
      ggplot2::geom_step(ggplot2::aes(
        x = .data$position / 1e6, y = .data$pct, colour = .data$group
      )) +
      ggplot2::facet_grid(~ chromosome, scales = "free_x", space = "free_x") +
      ggplot2::labs(
        x = "Position (Mb)", y = "Individuals in a run (%)",
        colour = "Group", title = run_word
      )
  } else if (type == "froh") {
    if (froh_by == "individual") {
      # F_ROH over the whole genome, for every individual
      froh <- .roh_total_length(runs, individuals, "sum")
      froh$froh <- froh$length_bp / sum(chromosome_length$length_bp)
      p <- ggplot2::ggplot(froh, ggplot2::aes(
        x = .data$group, y = .data$froh, colour = .data$group
      )) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.2, height = 0) +
        ggplot2::labs(
          x = "Group", y = "F_ROH", colour = "Group",
          title = "Inbreeding from runs of homozygosity (F_ROH)"
        )
    } else {
      # mean F_ROH of each group on each chromosome (individuals with no
      # run on a chromosome count as 0)
      n_indiv <- table(individuals$group)
      froh <- runs %>%
        dplyr::group_by(.data$group, .data$chromosome) %>%
        dplyr::summarise(
          length_bp = sum(.data$length_bp), .groups = "drop"
        ) %>%
        dplyr::left_join(chromosome_length,
          by = "chromosome", suffix = c("", "_chromosome")
        )
      froh$froh <- froh$length_bp / froh$length_bp_chromosome /
        as.numeric(n_indiv[as.character(froh$group)])
      p <- ggplot2::ggplot(froh, ggplot2::aes(
        x = .data$chromosome, y = .data$froh, fill = .data$group
      )) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::labs(
          x = "Chromosome", y = "Mean F_ROH", fill = "Group",
          title = "Inbreeding from runs of homozygosity (F_ROH) by chromosome"
        )
    }
  } else if (type == "violin") {
    totals <- .roh_total_length(runs, individuals, method)
    p <- ggplot2::ggplot(totals, ggplot2::aes(
      x = .data$group, y = .data$length_bp / 1e6, fill = .data$group
    )) +
      ggplot2::geom_violin() +
      ggplot2::geom_jitter(width = 0.1, height = 0, size = 0.5) +
      ggplot2::labs(
        x = "Group",
        y = paste(
          if (method == "sum") "Total" else "Mean",
          "length of runs per individual (Mb)"
        ),
        fill = "Group", title = run_word
      )
  } else if (type == "length_class") {
    if (!is.numeric(classes) || length(classes) < 1 || is.unsorted(classes)) {
      stop("classes should be increasing lower bounds of the classes, in Mb")
    }
    labels <- c(
      paste0(utils::head(classes, -1), "-", classes[-1]),
      paste0(">", classes[length(classes)])
    )
    runs$class <- cut(runs$length_bp / 1e6,
      breaks = c(classes, Inf), labels = labels, right = FALSE
    )
    pct <- runs %>%
      dplyr::filter(!is.na(.data$class)) %>%
      dplyr::count(.data$group, .data$class) %>%
      dplyr::group_by(.data$group) %>%
      dplyr::mutate(pct = 100 * .data$n / sum(.data$n)) %>%
      dplyr::ungroup()
    p <- ggplot2::ggplot(pct, ggplot2::aes(
      x = .data$class, y = .data$pct, fill = .data$group
    )) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(
        x = "Length of runs (Mb)", y = "Runs (%)", fill = "Group",
        title = run_word
      )
  }
  p + ggplot2::theme_minimal()
}

# The runs, with a `group` column (a factor): the label of the group of each
# run if the runs were found on a grouped gen_tibble (the grouping variables
# are the columns before `id`), otherwise "all"
.roh_runs_with_group <- function(object) {
  object <- tibble::as_tibble(object)
  object$group <- .roh_group_labels(object)
  object
}

# All the individuals (from the `individuals` attribute, or else those with
# runs), with their group as in .roh_runs_with_group(), ordered by group
.roh_individuals <- function(object) {
  individuals <- attr(object, "individuals")
  if (is.null(individuals)) {
    individuals <- tibble::as_tibble(object)[
      ,
      seq_len(match("id", names(object)))
    ]
    individuals <- individuals[!duplicated(individuals$id), ]
  }
  individuals$group <- .roh_group_labels(individuals)
  individuals$id <- as.character(individuals$id)
  individuals[order(individuals$group), c("id", "group")]
}

# the label of the group of each row (columns before `id`), as a factor with
# the groups in order
.roh_group_labels <- function(x) {
  group_cols <- names(x)[seq_len(match("id", names(x)) - 1)]
  if (length(group_cols) == 0) {
    return(factor(rep("all", nrow(x))))
  }
  if (length(group_cols) == 1 && is.factor(x[[group_cols]])) {
    return(droplevels(x[[group_cols]]))
  }
  keys <- dplyr::arrange(dplyr::distinct(x[group_cols]),
    dplyr::across(dplyr::everything())
  )
  factor(.group_labels(x[group_cols]), levels = .group_labels(keys))
}

# The total (or mean) length of the runs of every individual, with 0 for
# individuals with no runs
.roh_total_length <- function(runs, individuals, method) {
  per_indiv <- runs %>%
    dplyr::mutate(id = as.character(.data$id)) %>%
    dplyr::group_by(.data$id) %>%
    dplyr::summarise(length_bp = if (method == "sum") {
      sum(.data$length_bp)
    } else {
      mean(.data$length_bp)
    })
  res <- individuals
  res$length_bp <- per_indiv$length_bp[match(res$id, per_indiv$id)]
  res$length_bp[is.na(res$length_bp)] <- 0
  res
}

# The percentage of the individuals of each group that are in a run, along
# each chromosome: a step function that changes at the start and after the
# end of each run
.roh_coverage <- function(runs, individuals) {
  n_indiv <- table(individuals$group)
  events <- dplyr::bind_rows(
    tibble::tibble(
      group = runs$group, chromosome = runs$chromosome,
      position = runs$start, change = 1
    ),
    tibble::tibble(
      group = runs$group, chromosome = runs$chromosome,
      position = runs$end, change = -1
    )
  )
  events %>%
    dplyr::arrange(.data$group, .data$chromosome, .data$position,
      dplyr::desc(.data$change)
    ) %>%
    dplyr::group_by(.data$group, .data$chromosome) %>%
    dplyr::mutate(
      n = cumsum(.data$change),
      pct = 100 * .data$n / as.numeric(n_indiv[as.character(.data$group)])
    ) %>%
    dplyr::ungroup()
}
