#' Create a Quality Control report for loci
#'
#' Return QC information to assess loci (MAF, missingness and HWE test).
#'
#' @param .x a [`gen_tibble`] object.
#' @param ... currently unused
#' the HWE test.
#' @returns a tibble with 3 elements: maf, missingness and hwe_p
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
  stopifnot_diploid(.x$genotypes)
  message(paste(
    "This gen_tibble is not grouped. For Hardy-Weinberg equilibrium,",
    "`qc_report_loci()` will assume individuals are part of the same",
    "population and HWE test p-values will be calculated across all",
    "individuals. If you wish to calculate HWE p-values within",
    "populations or groups, please use`group_by()` before calling",
    "`qc_report_loci()`."
  ))
  qc_report <- .x %>%
    reframe(
      snp_id = loci_names(.x),
      maf = loci_maf(.data$genotypes),
      missingness = loci_missingness(.data$genotypes),
      hwe_p = loci_hwe(.data$genotypes)
    )
  class(qc_report) <- c("qc_report_loci", class(qc_report))
  qc_report
}

#' @export
#' @rdname qc_report_loci
qc_report_loci.grouped_df <- function(.x, ...) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x$genotypes)

  # Find number of groups
  pops <- .x %>%
    select(dplyr::group_vars(.x)) %>%
    dplyr::pull(1)

  # Take grouped variable pops and find how many unique groups
  pops <- length(unique(pops))

  # Calculate hwe across groups
  hwe_res <- loci_hwe(.x, type = "matrix")
  hwe_res <- as.data.frame(hwe_res)
  hwe_res$p_corrected <- apply(hwe_res, 1, function(row) min(row) * pops)

  qc_report <- .x %>%
    ungroup() %>%
    reframe(
      snp_id = loci_names(.x),
      maf = loci_maf(.data$genotypes),
      missingness = loci_missingness(.data$genotypes),
      hwe_p = hwe_res$p_corrected
    )
  class(qc_report) <- c("qc_report_loci", class(qc_report))
  qc_report
}

#' Autoplots for `qc_report_loci` objects
#'
#' For `qc_report_loci`, the following types of plots are available:
#' - `overview`: an UpSet plot, giving counts of snps over the threshold for
#' missingness, minor allele frequency, and Hardy-Weinberg equilibrium P-value,
#' and visualising the interaction between these
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
    hwe_p = 0.01, # SUGGESTION should this be hwe_p_threshold for consistency?
    ...) {
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

# TODO BUG autoplot


autoplot_l_qc_all <- function(
    object,
    maf_threshold,
    miss_threshold,
    hwe_p_low_thresh,
    hwe_p_vertical_line,
    ...) {
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
  return(combined_plots)
}


autoplot_l_qc_overview <- function(
    object,
    maf_threshold,
    miss_threshold,
    hwe_p_low_thresh,
    ...) {
  if (any(is.na(object))) {
    message(paste(
      "One or more loci are missing for every individual.",
      "These will be removed from the QC report plot."
    ))
    # remove NA's from object
    object <- object[!is.na(object$maf), ]
  }

  qc_report <- object

  qc_hwe <- qc_report[qc_report$hwe_p >= hwe_p_low_thresh, ]
  qc_maf <- qc_report[qc_report$maf >= maf_threshold, ]


  maf_pass <- c(qc_maf$snp_id)
  hwe_pass <- c(qc_hwe$snp_id)

  qc_missing <- qc_report[qc_report$missingness >= miss_threshold, ]
  missing_pass <- c(qc_missing$snp_id)

  pass_list <- list(MAF = maf_pass, HWE = hwe_pass, Missing = missing_pass)

  unique_markers <- unique(unlist(pass_list))
  pass_counts <- UpSetR::fromList(pass_list)
  rownames(pass_counts) <- unique_markers

  final_plot_overview <- UpSetR::upset(
    pass_counts,
    order.by = "freq",
    main.bar.color = "#66C2A5",
    matrix.color = "#66C2A5",
    sets.bar.color = "#FC8D62"
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
    ...) {
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
