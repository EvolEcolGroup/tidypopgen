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
qc_report_loci <- function(.x, ...) {
  UseMethod("qc_report_loci", .x)
}


#' @export
#' @rdname qc_report_loci
qc_report_loci.tbl_df <- function(.x, ...) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x$genotypes)
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
  hwe_res <- .x %>% group_map(.f = ~ loci_hwe(.x$genotypes))
  hwe_res <- do.call("cbind", hwe_res)
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
#' snps with low minor allele freqency
#' - `missing high maf`:a histogram of the proportion of missing data for
#' snps with high minor allele freqency
#' - `maf`: a histogram of minor allele frequency
#' - `hwe`: a histogram of HWE exact test p-values
#' - `significant hwe`: a histogram of significant HWE exact test p-values
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are
#' not customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' @param object an object of class `qc_report_loci`
#' @param type the type of plot (one of `overview`, `all`, `missing`,
#' `missing low maf`, `missing high maf`, `maf`, `hwe`, and `significant hwe`)
#' @param maf_threshold a threshold for the accepted rate of minor allele
#' frequency of loci
#' @param miss_threshold a threshold for the accepted rate of missingness per
#' loci
#' @param hwe_p a threshold of significance for Hardy-Weinberg exact p-values
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @export
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
  maf_threshold = NULL,
  miss_threshold = NULL,
  hwe_p = NULL,
  ...
) {
  type <- match.arg(type)

  rlang::check_dots_empty()

  maf_threshold <- if (is.null(maf_threshold)) {
    0.05
  } else {
    maf_threshold
  }

  miss_threshold <- if (is.null(miss_threshold)) {
    0.01
  } else {
    miss_threshold
  }

  hwe_p <- if (is.null(hwe_p)) {
    0.01
  } else {
    hwe_p
  }

  logp <- -log10(hwe_p)

  if (type == "overview") {
    final_plot <- autoplot_l_qc_overview(
      object,
      maf_threshold = maf_threshold,
      miss_threshold = miss_threshold,
      hwe_p = hwe_p
    )
  } else if (type == "all") {
    final_plot <-
      autoplot_l_qc_all(
        object,
        maf_threshold = maf_threshold,
        miss_threshold = miss_threshold,
        hwe_p = hwe_p,
        logp = logp
      )
  } else if (type == "missing") {
    final_plot <- autoplot_l_qc_missing(object, miss_threshold = miss_threshold)
  } else if (type == "missing low maf") {
    final_plot <-
      autoplot_l_qc_missing_low_maf(
        object,
        maf_threshold = maf_threshold,
        miss_threshold = miss_threshold
      )
  } else if (type == "missing high maf") {
    final_plot <-
      autoplot_l_qc_missing_high_maf(
        object,
        maf_threshold = maf_threshold,
        miss_threshold = miss_threshold
      )
  } else if (type == "maf") {
    final_plot <- autoplot_l_qc_maf(object, maf_threshold = maf_threshold)
  } else if (type == "hwe") {
    final_plot <- autoplot_l_qc_hwe(object, logp = logp)
  } else if (type == "significant hwe") {
    final_plot <- autoplot_l_qc_sig_hwe(object, hwe_p = hwe_p, logp = logp)
  } else {
    stop(paste(
      "Invalid type argument. Please choose from 'overview',",
      "'all','maf','hwe','significant hwe'"
    ))
  }

  return(final_plot)
}


autoplot_l_qc_all <- function(
  object,
  maf_threshold = maf_threshold,
  miss_threshold = miss_threshold,
  hwe_p = hwe_p,
  logp = logp,
  ...
) {
  qc_report <- object

  # Missingness (according to MAF thresholds)
  qc_highmaf <- subset(qc_report, qc_report$maf > maf_threshold)
  qc_lowmaf <- subset(qc_report, qc_report$maf <= maf_threshold)

  p_highmaf <- ggplot2::ggplot(
    qc_highmaf,
    ggplot2::aes(x = .data$missingness)
  ) +
    ggplot2::geom_histogram(
      position = "dodge",
      binwidth = 0.005,
      fill = "#66C2A5"
    ) +
    ggplot2::labs(
      x = "Proportion of missing data",
      y = "Number of SNPs",
      title = paste("SNPs with MAF > ", maf_threshold)
    ) +
    ggplot2::geom_vline(xintercept = miss_threshold, lty = 2, col = "red") +
    ggplot2::scale_color_brewer(palette = "Dark2")

  p_lowmaf <- ggplot2::ggplot(qc_lowmaf, ggplot2::aes(x = .data$missingness)) +
    ggplot2::geom_histogram(
      position = "dodge",
      binwidth = 0.005,
      fill = "#66C2A5"
    ) +
    ggplot2::labs(
      x = "Proportion of missing data",
      y = "Number of SNPs",
      title = paste("SNPs with MAF < ", maf_threshold)
    ) +
    ggplot2::geom_vline(xintercept = miss_threshold, lty = 2, col = "red") +
    ggplot2::scale_color_brewer(palette = "Dark2")

  mafmiss <- patchwork::wrap_plots(p_highmaf + p_lowmaf)

  # Hardy weinberg exact test p-val distribution
  qc_report$hwe_p_log <- -log10(qc_report$hwe_p)
  qc_lowhwe <- qc_report[qc_report$hwe_p < hwe_p, ]

  hwe_all <- ggplot2::ggplot(qc_report, ggplot2::aes(x = .data$hwe_p_log)) +
    ggplot2::geom_histogram(binwidth = 0.5, fill = "#66C2A5") +
    ggplot2::labs(
      x = expression("-log"[10] * " of HWE exact p-value"),
      y = "Number of SNPs",
      title = "Hardy-Weinberg exact test"
    ) +
    ggplot2::geom_vline(xintercept = logp, lty = 2, col = "red")

  hwe_low <- ggplot2::ggplot(qc_lowhwe, ggplot2::aes(x = .data$hwe_p_log)) +
    ggplot2::geom_histogram(binwidth = 0.5, fill = "#66C2A5") +
    ggplot2::labs(
      x = expression("-log"[10] * " of HWE exact p-value"),
      y = "Number of SNPs",
      title = paste(
        "HWE exact p-value <",
        hwe_p
      )
    ) +
    ggplot2::geom_vline(xintercept = logp, lty = 2, col = "red")

  hwes <- patchwork::wrap_plots(hwe_all, hwe_low)

  final_plot_all <- mafmiss / hwes
  return(final_plot_all)
}


autoplot_l_qc_overview <- function(
  object,
  maf_threshold = maf_threshold,
  miss_threshold = miss_threshold,
  hwe_p = hwe_p,
  ...
) {
  qc_report <- object

  qc_hwe <- qc_report[qc_report$hwe_p >= hwe_p, ]
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


autoplot_l_qc_maf <- function(object, maf_threshold = maf_threshold, ...) {
  qc_report <- object
  # Minor allele frequency distribution
  maf <- ggplot2::ggplot(qc_report, ggplot2::aes(x = .data$maf)) +
    ggplot2::geom_histogram(binwidth = 0.01, fill = "#66C2A5") +
    ggplot2::labs(
      x = "Minor allele frequency",
      y = "Number of SNPs",
      title = "Minor allele frequency distribution"
    ) +
    ggplot2::geom_vline(xintercept = maf_threshold, lty = 2, col = "red")
  return(maf)
}

autoplot_l_qc_hwe <- function(object, logp, ...) {
  qc_report <- object
  # Hardy weinberg exact test p-val distribution
  qc_report$hwe_p_log <- -log10(qc_report$hwe_p)
  hwe_all <- ggplot2::ggplot(qc_report, ggplot2::aes(x = .data$hwe_p_log)) +
    ggplot2::geom_histogram(binwidth = 0.5, fill = "#66C2A5") +
    ggplot2::labs(
      x = expression("-log"[10] * " of HWE exact p-value"),
      y = "Number of SNPs",
      title = "Hardy-Weinberg exact test"
    ) +
    ggplot2::geom_vline(xintercept = logp, lty = 2, col = "red")
  return(hwe_all) # nolint
}

autoplot_l_qc_sig_hwe <- function(object, hwe_p = hwe_p, logp, ...) {
  qc_report <- object
  # Hardy weinberg exact test p-val distribution
  qc_report$hwe_p_log <- -log10(qc_report$hwe_p)
  qc_lowhwe <- qc_report[qc_report$hwe_p < hwe_p, ]
  hwe_low <- ggplot2::ggplot(qc_lowhwe, ggplot2::aes(x = .data$hwe_p_log)) +
    ggplot2::geom_histogram(binwidth = 0.5, fill = "#66C2A5") +
    ggplot2::labs(
      x = expression("-log"[10] * " of HWE exact p-value"),
      y = "Number of SNPs",
      title = paste(
        "HWE exact p-value <",
        hwe_p
      )
    ) +
    ggplot2::geom_vline(xintercept = logp, lty = 2, col = "red")
  return(hwe_low)
}

autoplot_l_qc_missing <- function(
  object,
  miss_threshold = miss_threshold,
  ...
) {
  qc_report <- object
  missing <- ggplot2::ggplot(qc_report, ggplot2::aes(x = .data$missingness)) +
    ggplot2::geom_histogram(
      position = "dodge",
      binwidth = 0.005,
      fill = "#66C2A5"
    ) +
    ggplot2::labs(
      x = "Proportion of missing data",
      y = "Number of SNPs",
      title = "SNP missingness"
    ) +
    ggplot2::geom_vline(xintercept = miss_threshold, lty = 2, col = "red") +
    ggplot2::scale_color_brewer(palette = "Dark2")
  return(missing)
}


autoplot_l_qc_missing_low_maf <- function(
  object,
  maf_threshold = maf_threshold,
  miss_threshold = miss_threshold,
  ...
) {
  qc_report <- object

  qc_lowmaf <- subset(qc_report, qc_report$maf <= maf_threshold)

  p_highmaf <- ggplot2::ggplot(qc_lowmaf, ggplot2::aes(x = .data$missingness)) + # nolint
    ggplot2::geom_histogram(
      position = "dodge",
      binwidth = 0.005,
      fill = "#66C2A5"
    ) +
    ggplot2::labs(
      x = "Proportion of missing data",
      y = "Number of SNPs",
      title = paste("SNPs with MAF <", maf_threshold)
    ) +
    ggplot2::geom_vline(xintercept = miss_threshold, lty = 2, col = "red") +
    ggplot2::scale_color_brewer(palette = "Dark2")
}

autoplot_l_qc_missing_high_maf <- function(
  object,
  maf_threshold = maf_threshold,
  miss_threshold = miss_threshold,
  ...
) {
  qc_report <- object

  qc_highmaf <- subset(qc_report, qc_report$maf > maf_threshold)

  p_highmaf <- ggplot2::ggplot(
    # nolint
    qc_highmaf,
    ggplot2::aes(x = .data$missingness)
  ) +
    ggplot2::geom_histogram(
      position = "dodge",
      binwidth = 0.005,
      fill = "#66C2A5"
    ) +
    ggplot2::labs(
      x = "Proportion of missing data",
      y = "Number of SNPs",
      title = paste("SNPs with MAF >", maf_threshold)
    ) +
    ggplot2::geom_vline(xintercept = miss_threshold, lty = 2, col = "red") +
    ggplot2::scale_color_brewer(palette = "Dark2")
}
