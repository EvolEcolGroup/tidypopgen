#' Create a Quality Control report for loci
#'
#' Return QC information to assess loci (MAF, missingness and HWE test). For
#' pseudohaploid data, HWE test is not calculated.
#'
#' @param .x a [`gen_tibble`] object.
#' @param ... currently unused
#' @returns either a tibble with 3 elements (maf, missingness and hwe_p). For
#'   pseudohaploid data, a tibble with 2 elements (maf and missingness).
#' @rdname qc_report_loci
#' @export
#' @examples
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' example_gt <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Get a QC report for the loci
#' example_gt %>% qc_report_loci()
#'
#' # Group by population to calculate HWE within populations
#' example_gt <- example_gt %>% group_by(population)
#' example_gt %>% qc_report_loci()
#'
qc_report_loci <- function(.x, ...) {
  UseMethod("qc_report_loci", .x)
}

#' @export
#' @rdname qc_report_loci
qc_report_loci.tbl_df <- function(.x, ...) {
  rlang::check_dots_empty()

  if (show_ploidy(.x) == -2) {
    qc_report <- qc_report_loci_pseudohaploid(.x)
    return(qc_report)
  }

  stopifnot_diploid(.x$genotypes)
  message(paste(
    "This gen_tibble is not grouped. For Hardy-Weinberg equilibrium,",
    "`qc_report_loci()` will assume individuals are part of the same",
    "population and HWE test p-values will be calculated across all",
    "individuals. If you wish to calculate HWE p-values within",
    "populations or groups, please use`group_by()` before calling",
    "`qc_report_loci()`."
  ))
  # compute maf, missingness and hwe (treating all individuals as a single
  # group)
  qc_report <- qc_report_loci_stats(.x, hwe = TRUE)
  class(qc_report) <- c("qc_report_loci", class(qc_report))
  return(qc_report)
}

qc_report_loci_pseudohaploid <- function(.x, ...) {
  rlang::check_dots_empty()
  # maf and missingness are computed across all individuals, without HWE
  qc_report <- qc_report_loci_stats(ungroup(.x), hwe = FALSE)
  class(qc_report) <- c("qc_report_loci", class(qc_report))
  return(qc_report)
}

#' Compute the statistics for qc_report_loci
#'
#' Internal function that computes the minor allele frequency and missingness
#' of each locus across all individuals and, if `hwe = TRUE`, the HWE exact
#' test p-value (mid-p) within each group, returning the smallest p-value
#' across groups (a single group for an ungrouped `gen_tibble`). When
#' genotypes are coded as 0, 1, 2 or missing (raw or imputed calls), all
#' statistics are computed in a single pass over the genotypes by
#' `loci_qc_cpp()`, rather than reading the genotypes once for each
#' statistic; otherwise (e.g. dosages), or for a single individual, it falls
#' back on [loci_maf()], [loci_missingness()] and [loci_hwe()].
#'
#' @param .x a [`gen_tibble`] object, possibly grouped.
#' @param hwe boolean, whether the HWE test should be computed.
#' @returns a tibble with columns `snp_id`, `maf`, `missingness` and (if
#'   `hwe = TRUE`) `hwe_p`, the smallest p-value across groups.
#' @keywords internal
#' @noRd
qc_report_loci_stats <- function(.x, hwe) {
  geno_fbm <- .gt_get_fbm(.x)
  # group of each individual (all in one group if not grouped)
  if (inherits(.x, "grouped_df")) {
    group_ids <- dplyr::group_indices(.x)
    n_groups <- dplyr::n_groups(.x)
  } else {
    group_ids <- rep(1L, nrow(.x))
    n_groups <- 1L
  }
  code256 <- geno_fbm$code256
  if (all(is.na(code256) | code256 %in% c(0, 1, 2)) && nrow(.x) > 1) {
    # a single pass over the genotypes (openMP is used for parallelisation)
    qc_stats <- loci_qc_cpp(
      BM = geno_fbm,
      rowInd = .gt_fbm_rows(.x),
      colInd = .gt_fbm_cols(.x),
      groupIds = group_ids - 1L,
      ngroups = n_groups,
      ploidy = indiv_ploidy(.x),
      hwe = hwe,
      midp = TRUE,
      ncores = bigstatsr::nb_cores(),
      block_size = bigstatsr::block_size(nrow(.x), 1)
    )
  } else {
    # compute each statistic separately
    qc_stats <- list(
      maf = loci_maf(.x$genotypes),
      missingness = loci_missingness(.x$genotypes)
    )
    if (hwe) {
      if (inherits(.x, "grouped_df")) {
        # smallest p-value across groups (one vector per group)
        qc_stats$hwe_p <- do.call(pmin, loci_hwe(.x, type = "list"))
      } else {
        qc_stats$hwe_p <- loci_hwe(.x$genotypes)
      }
    }
  }
  qc_report <- tibble::tibble(
    snp_id = loci_names(.x),
    maf = qc_stats$maf,
    missingness = qc_stats$missingness
  )
  if (hwe) {
    qc_report$hwe_p <- qc_stats$hwe_p
  }
  qc_report
}

#' @export
#' @rdname qc_report_loci
qc_report_loci.grouped_df <- function(.x, ...) {
  rlang::check_dots_empty()

  if (show_ploidy(.x) == -2) {
    qc_report <- qc_report_loci_pseudohaploid(.x)
    return(qc_report)
  }

  stopifnot_diploid(.x$genotypes)

  # Find the number of groups tested for HWE (empty groups, which only exist
  # if the user asked for them, are not tested)
  pops <- sum(lengths(dplyr::group_rows(.x)) > 0)

  # Calculate maf and missingness across all individuals, and the smallest
  # hwe p-value across groups (in a single pass over the genotypes)
  qc_report <- qc_report_loci_stats(.x, hwe = TRUE)
  # Bonferroni correction of the smallest p-value, capped at 1
  qc_report$hwe_p <- pmin(1, qc_report$hwe_p * pops)
  class(qc_report) <- c("qc_report_loci", class(qc_report))
  return(qc_report)
}

#' Autoplots for `qc_report_loci` objects
#'
#' For `qc_report_loci`, the following types of plots are available:
#' - `overview`: an UpSet plot, giving counts of snps over the threshold for
#' missingness, minor allele frequency, and Hardy-Weinberg equilibrium P-value,
#' and visualising the interaction between these.
#' - `all`: a four panel plot, containing `missing high maf`, `missing low maf`,
#' `hwe`, and `significant hwe` plots
#' - `missing`: a histogram of proportion of missing data
#' - `missing low maf`: a histogram of the proportion of missing data for
#' snps with low minor allele frequency
#' - `missing high maf`:a histogram of the proportion of missing data for
#' snps with high minor allele frequency
#' - `maf`: a histogram of minor allele frequency
#' - `hwe`: a histogram of HWE exact test p-values
#' - `significant hwe`: a histogram of significant HWE exact test p-values
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are not
#' customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' @param object an object of class `qc_report_loci`
#' @param type the type of plot (one of `overview`, `all`, `missing`,
#' `missing low maf`, `missing high maf`, `maf`, `hwe`, and `significant hwe`)
#' @param maf_threshold default 0.05, a threshold for the accepted rate of minor
#'   allele frequency of loci
#' @param miss_threshold default 0.01, a threshold for the accepted rate of
#'   missingness per loci
#' @param hwe_p default 0.01, a threshold of significance for Hardy-Weinberg
#'   exact p-values
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @export
#' @examples
#' # Create a gen_tibble
#' bed_file <-
#'   system.file("extdata", "related", "families.bed", package = "tidypopgen")
#' example_gt <- gen_tibble(bed_file,
#'   backingfile = tempfile("families"),
#'   quiet = TRUE,
#'   valid_alleles = c("1", "2")
#' )
#'
#' loci_report <- example_gt %>% qc_report_loci()
#'
#' # Plot the QC report overview
#' autoplot(loci_report, type = "overview")
#'
#' # Plot the QC report all
#' autoplot(loci_report, type = "all")
#'
#' # Plot missing data
#' autoplot(loci_report, type = "missing")
#'
#' # Plot missing with low maf
#' autoplot(loci_report, type = "missing low maf", maf_threshold = 0.05)
#'
#' # Plot missing with high maf
#' autoplot(loci_report, type = "missing high maf", maf_threshold = 0.05)
#'
#' # Plot maf
#' autoplot(loci_report, type = "maf", maf_threshold = 0.05)
#'
#' # Plot hwe
#' autoplot(loci_report, type = "hwe", hwe_p = 0.01)
#'
#' # Plot significant hwe
#' autoplot(loci_report, type = "significant hwe", hwe_p = 0.01)
#'
autoplot.qc_report_loci <- function(
  object,
  type = c(
    "overview",
    "all",
    "missing",
    "missing low maf",
    "missing high maf",
    "maf",
    "hwe",
    "significant hwe"
  ),
  maf_threshold = 0.05,
  miss_threshold = 0.01,
  hwe_p = 0.01,
  ...
) {
  type <- match.arg(type)

  rlang::check_dots_empty()

  if (type == "overview") {
    report_plot <- autoplot_l_qc_overview(
      object,
      maf_threshold = maf_threshold,
      miss_threshold = miss_threshold,
      hwe_p_low_thresh = hwe_p
    )
  } else if (type == "all") {
    report_plot <-
      autoplot_l_qc_all(
        object,
        maf_threshold = maf_threshold,
        miss_threshold = miss_threshold,
        hwe_p_low_thresh = hwe_p,
        hwe_p_vertical_line = hwe_p
      )
  } else if (type == "missing") {
    report_plot <- autoplot_l_qc_missing(object,
      miss_threshold = miss_threshold
    )
  } else if (type == "missing low maf") {
    report_plot <-
      autoplot_l_qc_missing(
        object,
        maf_low_thresh = maf_threshold,
        miss_threshold = miss_threshold
      )
  } else if (type == "missing high maf") {
    report_plot <-
      autoplot_l_qc_missing(
        object,
        maf_high_thresh = maf_threshold,
        miss_threshold = miss_threshold
      )
  } else if (type == "maf") {
    report_plot <- autoplot_l_qc_maf(object, maf_threshold = maf_threshold)
  } else if (type == "hwe") {
    report_plot <- autoplot_l_qc_hwe(object, hwe_p_vertical_line = hwe_p)
  } else if (type == "significant hwe") {
    report_plot <- autoplot_l_qc_hwe(object,
      hwe_p_low_thresh = hwe_p,
      hwe_p_vertical_line = hwe_p
    )
  }

  return(report_plot)
}

autoplot_l_qc_all <- function(
  object,
  maf_threshold,
  miss_threshold,
  hwe_p_low_thresh,
  hwe_p_vertical_line,
  ...
) {
  # Missingness (according to MAF thresholds)
  miss_high_maf_plot <- autoplot_l_qc_missing(
    object,
    miss_threshold = miss_threshold,
    maf_high_thresh = maf_threshold
  )
  miss_low_maf_plot <- autoplot_l_qc_missing(
    object,
    miss_threshold = miss_threshold,
    maf_low_thresh = maf_threshold
  )

  miss_maf_plots <- patchwork::wrap_plots(
    miss_high_maf_plot +
      miss_low_maf_plot
  )


  if ("hwe_p" %in% colnames(object)) {
    # Hardy Weinberg exact test p-val distribution
    hwe_all_plot <- autoplot_l_qc_hwe(object,
      hwe_p_vertical_line = hwe_p_vertical_line
    )
    hwe_low_plot <- autoplot_l_qc_hwe(
      object = object,
      hwe_p_vertical_line = hwe_p_vertical_line,
      hwe_p_low_thresh = hwe_p_low_thresh
    )

    hwe_plots <- patchwork::wrap_plots(hwe_all_plot, hwe_low_plot)

    combined_plots <- miss_maf_plots / hwe_plots
  } else {
    combined_plots <- miss_maf_plots
  }
  return(combined_plots)
}


autoplot_l_qc_overview <- function(
  object,
  maf_threshold,
  miss_threshold,
  hwe_p_low_thresh,
  ...
) {
  # check that the ellipses are empty
  if (length(list(...)) > 0) {
    stop("no additional arguments should be provided for this plot")
  }
  if (any(is.na(object))) {
    message(paste(
      "One or more loci are missing for every individual.",
      "These will be removed from the QC report plot."
    ))
    # remove NA's from object
    object <- object[!is.na(object$maf), ]
  }

  qc_report <- data.frame(
    Missing = object$missingness < miss_threshold,
    MAF = object$maf > maf_threshold
  )
  # if including HWE
  if ("hwe_p" %in% colnames(object)) {
    qc_report$HWE <- object$hwe_p >= hwe_p_low_thresh
  }
  final_plot_overview <- upset_plot(qc_report,
    bar_colour = "#66C2A5",
    dot_colour = "#66C2A5",
    set_bar_colour = "#FC8D62"
  )
  return(final_plot_overview)
}


autoplot_l_qc_maf <- function(object, maf_threshold, ...) {
  # Minor allele frequency distribution
  maf <- ggplot2::ggplot(object, ggplot2::aes(x = .data$maf)) +
    ggplot2::geom_histogram(binwidth = 0.01, fill = "#66C2A5") +
    ggplot2::labs(
      x = "Minor allele frequency",
      y = "Number of SNPs",
      title = "Minor allele frequency distribution"
    ) +
    ggplot2::geom_vline(xintercept = maf_threshold, lty = 2, col = "red")
  return(maf)
}

#' Internal plotting function for hwe
#'
#' @param object an object of class `qc_report_loci`
#' @param hwe_p_vertical_line the p-value threshold for the vertical line
#' @param hwe_p_low_thresh the p-value threshold to subset the loci
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @keywords internal
#' @noRd
autoplot_l_qc_hwe <- function(object,
                              hwe_p_vertical_line,
                              hwe_p_low_thresh = NULL,
                              ...) {
  # check ellipses are empty
  rlang::check_dots_empty()

  if (!"hwe_p" %in% colnames(object)) {
    stop(
      "Hardy-Weinberg equilibrium test not available for",
      "pseudohaploid data."
    )
  }

  hwe_p_vertical_line <- -log10(hwe_p_vertical_line)

  object$hwe_p_log <- -log10(object$hwe_p)
  # subset if necessary
  if (!is.null(hwe_p_low_thresh)) {
    object <- object[object$hwe_p < hwe_p_low_thresh, ]
    plot_title <- paste("HWE exact p-value <", hwe_p_low_thresh)
  } else {
    plot_title <- "Hardy-Weinberg exact test"
  }
  # Hardy Weinberg exact test p-val distribution
  hwe_plot <- ggplot2::ggplot(object, ggplot2::aes(x = .data$hwe_p_log)) +
    ggplot2::geom_histogram(binwidth = 0.5, fill = "#66C2A5") +
    ggplot2::labs(
      x = expression("-log"[10] * " of HWE exact p-value"),
      y = "Number of SNPs",
      title = plot_title
    ) +
    ggplot2::geom_vline(xintercept = hwe_p_vertical_line, lty = 2, col = "red")
  return(hwe_plot)
}

#' Internal missingness plot function
#'
#' @param object an object of class `qc_report_loci`
#' @param miss_threshold a threshold for the accepted rate of missingness per
#' loci
#' @param maf_low_thresh a threshold for the accepted rate of minor allele
#' frequency of loci
#' @param maf_high_thresh a threshold for the accepted rate of minor allele
#' frequency of loci
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @keywords internal
#' @noRd
autoplot_l_qc_missing <- function(
  object,
  miss_threshold,
  maf_low_thresh = NULL,
  maf_high_thresh = NULL,
  ...
) {
  # check ellipses are empty
  rlang::check_dots_empty()

  # subset if necessary
  if (!is.null(maf_low_thresh)) {
    object <- subset(object, object$maf <= maf_low_thresh)
    plot_title <- paste("SNPs with MAF <", maf_low_thresh)
  } else if (!is.null(maf_high_thresh)) {
    object <- subset(object, object$maf > maf_high_thresh)
    plot_title <- paste("SNPs with MAF >", maf_high_thresh)
  } else {
    plot_title <- "SNP missingness"
  }

  missing_plot <- ggplot2::ggplot(object, ggplot2::aes(x = .data$missingness)) +
    ggplot2::geom_histogram(
      position = "dodge",
      binwidth = 0.005,
      fill = "#66C2A5"
    ) +
    ggplot2::labs(
      x = "Proportion of missing data",
      y = "Number of SNPs",
      title = plot_title
    ) +
    ggplot2::geom_vline(xintercept = miss_threshold, lty = 2, col = "red") +
    ggplot2::scale_color_brewer(palette = "Dark2")
  return(missing_plot)
}
