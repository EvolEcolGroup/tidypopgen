#' Plot site frequency spectra
#'
#' Plot the site frequency spectra computed by [pop_sfs()] (bar plots of the
#' number of sites by allele count, one panel per population),
#' [pairwise_pop_sfs()] (heatmaps of the joint spectra, one panel per pair of
#' populations), or a single spectrum returned with `type = "array"` (a bar
#' plot for one population, a heatmap for two). The output of
#' [pop_sfs_preview()] is plotted as the number of segregating sites for each
#' projection, one panel per population, with the projection chosen by
#' default marked. Spectra of more than two populations (from
#' [nwise_pop_sfs()]) cannot be plotted; use [pairwise_pop_sfs()] to plot the
#' joint spectra of pairs of populations.
#'
#' The masked cells of the spectra (the corners, i.e. sites monomorphic in
#' the projected sample, and the cells that were folded out in a folded
#' spectrum; see [pop_sfs()]) are not shown by default. In heatmaps, the
#' numbers of sites are shown on a log scale by default, and empty cells
#' are left blank. As the projection spreads each locus over many cells, many
#' cells have a tiny fraction of a site: on a log scale, the colour scale
#' starts at `min_sites` (by default, a thousandth of the largest number of
#' sites in a cell), and cells with fewer sites are shown with its lightest
#' colour (as `vmin` in the plots of `dadi`).
#'
#' @param object the output of [pop_sfs()] or [pairwise_pop_sfs()] (with
#'   `type = "tidy"`, the default), an `sfs` object of one or two populations
#'   (returned with `type = "array"`), or the output of [pop_sfs_preview()].
#' @param masked boolean, whether to show the masked cells of the spectra.
#'   Defaults to `FALSE`.
#' @param log boolean, whether to show the numbers of sites on a log scale
#'   (the y axis of bar plots, the colour scale of heatmaps). If `NULL` (the
#'   default), a log scale is used for heatmaps only.
#' @param min_sites the lower limit of the log colour scale of heatmaps:
#'   cells with fewer (but more than 0) sites are shown with the lightest
#'   colour. If `NULL` (the default), a thousandth of the largest number of
#'   sites in a cell.
#' @param ... currently unused.
#' @returns a `ggplot2` object.
#' @seealso [pop_sfs()], [pairwise_pop_sfs()], [pop_sfs_preview()]
#' @examples
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # the 1D spectrum of each population
#' example_gt %>%
#'   pop_sfs(n_proj = 4) %>%
#'   autoplot()
#'
#' # the joint spectrum of each pair of populations
#' example_gt %>%
#'   pairwise_pop_sfs(n_proj = 4, folded = FALSE) %>%
#'   autoplot()
#'
#' # a single spectrum
#' sfs_2d <- example_gt %>%
#'   pairwise_pop_sfs(n_proj = 4, folded = FALSE, type = "array")
#' autoplot(sfs_2d$pop1.pop2)
#'
#' # the number of segregating sites for each projection
#' example_gt %>%
#'   pop_sfs_preview() %>%
#'   autoplot()
#' @name autoplot_sfs
#' @aliases autoplot_pop_sfs_tbl autoplot_pairwise_sfs_tbl
#' @aliases autoplot_sfs_preview_tbl
#' @importFrom rlang .data
#' @export
autoplot.sfs <- function(object, masked = FALSE, log = NULL,
                         min_sites = NULL, ...) {
  rlang::check_dots_empty()
  .sfs_check_flag(masked, "masked")
  pops <- attr(object, "pop_ids")
  folded <- attr(object, "folded")
  if (length(pops) == 1) {
    cells <- .sfs_cells(object, count_names = "allele_count")
    cells$panel <- factor(pops, levels = pops)
    .sfs_plot_1d(cells,
      folded = folded, masked = masked, log = isTRUE(log),
      facet = FALSE
    )
  } else if (length(pops) == 2) {
    cells <- .sfs_cells(object,
      count_names = c("allele_count_1", "allele_count_2")
    )
    cells$panel <- factor(paste(pops, collapse = "."))
    .sfs_plot_2d(cells,
      folded = folded, masked = masked, log = !isFALSE(log),
      min_sites = min_sites, facet = FALSE, axis_names = pops
    )
  } else {
    stop(
      "only spectra of one or two populations can be plotted; use ",
      "pairwise_pop_sfs() to plot the joint spectra of pairs of populations"
    )
  }
}

#' @rdname autoplot_sfs
#' @export
autoplot.pop_sfs_tbl <- function(object, masked = FALSE, log = NULL, ...) {
  rlang::check_dots_empty()
  .sfs_check_flag(masked, "masked")
  .sfs_check_columns(object, c("allele_count", "n_sites", "masked"), "pop_sfs")
  group_cols <- setdiff(names(object), c("allele_count", "n_sites", "masked"))
  object$panel <- .sfs_panel_labels(object, group_cols)
  .sfs_plot_1d(object,
    folded = isTRUE(attr(object, "folded")), masked = masked,
    log = isTRUE(log), facet = length(group_cols) > 0
  )
}

#' @rdname autoplot_sfs
#' @export
autoplot.pairwise_sfs_tbl <- function(object, masked = FALSE, log = NULL,
                                      min_sites = NULL, ...) {
  rlang::check_dots_empty()
  .sfs_check_flag(masked, "masked")
  count_cols <- c("allele_count_1", "allele_count_2", "n_sites", "masked")
  .sfs_check_columns(object, count_cols, "pairwise_pop_sfs")
  group_cols <- setdiff(names(object), count_cols)
  group_1 <- .sfs_paste_groups(object, grep("_1$", group_cols, value = TRUE))
  group_2 <- .sfs_paste_groups(object, grep("_2$", group_cols, value = TRUE))
  # label each pair as "pop1 (x) - pop2 (y)", keeping the order of the pairs
  object$panel <- .sfs_panel_factor(
    object[group_cols],
    paste0(group_1, " (x) - ", group_2, " (y)")
  )
  .sfs_plot_2d(object,
    folded = isTRUE(attr(object, "folded")), masked = masked,
    log = !isFALSE(log), min_sites = min_sites, facet = TRUE
  )
}

#' @rdname autoplot_sfs
#' @export
autoplot.sfs_preview_tbl <- function(object, ...) {
  rlang::check_dots_empty()
  preview_cols <- c("n_proj", "n_seg_sites", "n_loci", "best")
  .sfs_check_columns(object, preview_cols, "pop_sfs_preview")
  group_cols <- setdiff(names(object), preview_cols)
  object$panel <- .sfs_panel_labels(object, group_cols)
  best <- object[object$best, , drop = FALSE]
  best$label <- paste0("n_proj = ", best$n_proj)
  p <- ggplot2::ggplot(object, ggplot2::aes(
    x = .data$n_proj,
    y = .data$n_seg_sites
  )) +
    ggplot2::geom_line(colour = .sfs_colours$bar, linewidth = 0.75) +
    # the projection chosen by default, with a ring in the background colour
    ggplot2::geom_point(
      data = best, shape = 21, size = 3, stroke = 1,
      fill = .sfs_colours$bar, colour = "white"
    ) +
    ggplot2::geom_text(
      data = best, ggplot2::aes(label = .data$label),
      vjust = -1.2, size = 3, colour = .sfs_colours$ink
    ) +
    ggplot2::scale_x_continuous(breaks = .sfs_int_breaks) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.15))
    ) +
    ggplot2::labs(
      x = "Projection (number of allele copies)",
      y = "Number of segregating sites"
    ) +
    .sfs_theme()
  if (length(group_cols) > 0) {
    p <- p + ggplot2::facet_wrap(~ .data$panel, scales = "free")
  }
  p
}

# colours: a single blue ramp (light to dark) for the numbers of sites, and a
# dark grey for text
.sfs_colours <- list(
  bar = "#2a78d6",
  low = "#cde2fb",
  high = "#0d366b",
  ink = "#3d3d3d"
)

# a minimal theme with recessive grid lines
.sfs_theme <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# the argument setting the transformation of a scale, which was called trans
# before ggplot2 3.5.0
.sfs_transform_arg <- function(transform) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    list(transform = transform)
  } else {
    list(trans = transform)
  }
}

# numbers without scientific notation
.sfs_number_labels <- function(x) {
  format(x,
    scientific = FALSE, trim = TRUE, drop0trailing = TRUE,
    big.mark = ","
  )
}

# integer breaks for axes of allele counts
.sfs_int_breaks <- function(limits) {
  breaks <- pretty(limits)
  breaks[breaks == round(breaks)]
}

# check that a tidy spectrum has the columns needed to plot it
.sfs_check_columns <- function(object, columns, fun) {
  if (!all(columns %in% names(object))) {
    stop(
      "`object` should be the output of ", fun, "(), with columns ",
      paste(columns, collapse = ", ")
    )
  }
  invisible(TRUE)
}

# a label for each row, from the values of the group columns (pasted with
# "."), as a factor with the groups in the order in which they appear
.sfs_panel_labels <- function(object, group_cols) {
  if (length(group_cols) == 0) {
    return(factor(rep("all", nrow(object))))
  }
  .sfs_panel_factor(object[group_cols], .sfs_paste_groups(object, group_cols))
}

# the values of the group columns of each row, pasted with "."
.sfs_paste_groups <- function(object, group_cols) {
  do.call(paste, c(
    unname(lapply(object[group_cols], as.character)),
    sep = "."
  ))
}

# the panel of each row, as a factor with the panels in the order in which
# they appear: a panel is a combination of values of the group columns
# (`keys`), labelled with `labels`. Combinations with the same label (e.g.
# ("a.b", "c") and ("a", "b.c"), or NA and "NA") are kept in their own
# panels, with labels made unique
.sfs_panel_factor <- function(keys, labels) {
  # ids in the order in which the combinations appear
  panel_id <- vctrs::vec_group_id(keys)
  first_row <- match(seq_len(attr(panel_id, "n")), panel_id)
  panel_labels <- make.unique(labels[first_row], sep = " ")
  factor(panel_labels[panel_id], levels = panel_labels)
}

# bar plot of 1D spectra (cells with columns allele_count, n_sites, masked
# and panel)
.sfs_plot_1d <- function(cells, folded, masked, log, facet) {
  if (!masked) {
    cells <- cells[!cells$masked, , drop = FALSE]
  }
  if (log) {
    # empty cells cannot be shown on a log scale
    cells <- cells[cells$n_sites > 0, , drop = FALSE]
  }
  # e.g. with n_proj = 1, both cells are masked (and facet_wrap() fails
  # without data)
  if (nrow(cells) == 0) {
    stop(
      "there are no cells to plot (they are all masked",
      if (log) ", or empty", "); use `masked = TRUE` to show the masked cells"
    )
  }
  p <- ggplot2::ggplot(cells, ggplot2::aes(
    x = .data$allele_count,
    y = .data$n_sites
  )) +
    ggplot2::geom_col(fill = .sfs_colours$bar, width = 0.7) +
    ggplot2::scale_x_continuous(breaks = .sfs_int_breaks) +
    ggplot2::labs(
      x = if (folded) "Minor allele count" else "Derived allele count",
      y = "Number of sites"
    ) +
    .sfs_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  if (log) {
    p <- p + ggplot2::scale_y_log10()
  }
  if (facet) {
    p <- p + ggplot2::facet_wrap(~ .data$panel, scales = "free")
  } else {
    p <- p + ggplot2::ggtitle(as.character(cells$panel[1]))
  }
  p
}

# heatmap of 2D spectra (cells with columns allele_count_1, allele_count_2,
# n_sites, masked and panel); masked (unless masked = TRUE) and, on a log
# scale, empty cells are left blank, and cells with fewer than min_sites are
# shown with the lightest colour
.sfs_plot_2d <- function(cells, folded, masked, log, facet,
                         min_sites = NULL, axis_names = NULL) {
  if (!masked) {
    cells$n_sites[cells$masked] <- NA
  }
  limits <- NULL
  if (log) {
    cells$n_sites[cells$n_sites <= 0] <- NA
    if (any(!is.na(cells$n_sites))) {
      if (is.null(min_sites)) {
        min_sites <- max(cells$n_sites, na.rm = TRUE) / 1000
      }
      is_number <- is.numeric(min_sites) && length(min_sites) == 1 &&
        !is.na(min_sites) && min_sites > 0
      if (!is_number) {
        stop("`min_sites` should be a positive number")
      }
      cells$n_sites <- pmax(cells$n_sites, min_sites)
      limits <- c(min_sites, NA)
    }
  }
  allele <- if (folded) "Minor allele count" else "Derived allele count"
  if (is.null(axis_names)) {
    x_lab <- paste(allele, "(first population)")
    y_lab <- paste(allele, "(second population)")
  } else {
    x_lab <- paste0(allele, " (", axis_names[1], ")")
    y_lab <- paste0(allele, " (", axis_names[2], ")")
  }
  p <- ggplot2::ggplot(cells, ggplot2::aes(
    x = .data$allele_count_1,
    y = .data$allele_count_2,
    fill = .data$n_sites
  )) +
    # the white borders separate the cells
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    do.call(ggplot2::scale_fill_gradient, c(
      list(
        low = .sfs_colours$low,
        high = .sfs_colours$high,
        na.value = "transparent",
        limits = limits,
        labels = .sfs_number_labels,
        name = "Number\nof sites"
      ),
      .sfs_transform_arg(if (log) "log10" else "identity")
    )) +
    ggplot2::scale_x_continuous(breaks = .sfs_int_breaks, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = .sfs_int_breaks, expand = c(0, 0)) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    .sfs_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (facet) {
    p <- p + ggplot2::facet_wrap(~ .data$panel, scales = "free") +
      ggplot2::theme(aspect.ratio = 1)
  } else {
    p <- p + ggplot2::coord_equal() +
      ggplot2::ggtitle(as.character(cells$panel[1]))
  }
  p
}
