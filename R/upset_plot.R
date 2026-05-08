#' Draw an UpSet plot from a data frame of logical columns
#'
#' This is a minimalistic implementation of upset plots, as used to visualise
#' our qc_loci_reports.  It is not intended to be a general-purpose upset
#' plotting function, and it does not attempt to replicate all the features of
#' specialised packages. It produces a three-panel UpSet plot assembled with
#' Produces a three-panel UpSet plot assembled with
#' \code{\link[patchwork]{patchwork}}:
#' \itemize{
#'   \item \strong{Top panel} - bar chart of intersection sizes, labelled.
#'   \item \strong{Matrix panel} - dot-and-line membership matrix with
#'     alternating row backgrounds for readability.
#'   \item \strong{Left panel} - horizontal bar chart of per-set totals.
#' }
#'
#' @param df             A data frame with logical set-membership columns.
#' @param sets           Character vector of column names to use as sets.
#'   Defaults to \code{NULL}, in which case all \code{logical} columns are
#'   used automatically.  The order of \code{sets} controls the top-to-bottom
#'   row order in the matrix (first element appears at the top).
#' @param min_size       Integer.  Intersections with fewer than this many rows
#'   are dropped.  Default \code{1L}.
#' @param n_intersections Integer.  Maximum number of intersections to display
#'   (top-n by count).  Default \code{40L}.
#' @param bar_colour     Fill colour for the intersection size bars.
#'   Default \code{"#2166ac"}.
#' @param dot_colour     Colour for filled dots and connecting lines in the
#'   membership matrix.  Default \code{"#2166ac"}.
#' @param empty_colour   Colour for absent-set dots in the matrix.
#'   Default \code{"#d9d9d9"}.
#' @param set_bar_colour Fill colour for the set-size bars.
#'   Default \code{"#4dac26"}.
#' @param text_size      Base font size (pts) passed to all
#'   \code{\link[ggplot2]{theme}} calls.  Default \code{11}.
#'
#' @return A \code{\link[patchwork]{patchwork}} object.  Print it or pass it to
#'   \code{\link[ggplot2]{ggsave}}.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   A = sample(c(TRUE, FALSE), 100, replace = TRUE),
#'   B = sample(c(TRUE, FALSE), 100, replace = TRUE),
#'   C = sample(c(TRUE, FALSE), 100, replace = TRUE)
#' )
#' upset_plot(df)
#'
#' @export
upset_plot <- function(
    df,
    sets             = NULL,
    min_size         = 1L,
    n_intersections  = 40L,
    bar_colour       = "#2166ac",
    dot_colour       = "#2166ac",
    empty_colour     = "#d9d9d9",
    set_bar_colour   = "#4dac26",
    text_size        = 11
) {

  # ── 0. resolve set columns ─────────────────────────────────────────────────
  if (is.null(sets)) {
    sets <- names(df)[vapply(df, is.logical, logical(1L))]
    if (length(sets) == 0L)
      stop("No logical columns found in `df`. Supply `sets` explicitly.")
  }

  # Single authoritative factor level definition used by every panel.
  # Bottom of the y-axis = sets[1], top = sets[length(sets)].
  # scale_y_discrete() places level 1 at the bottom, so reversing `sets`
  # puts the first element at the top visually.
  set_levels <- rev(sets)

  # Named integer lookup: set name -> y-axis position (integer).
  # Used to compute segment endpoints independent of factor coercion order.
  set_pos <- stats::setNames(seq_along(set_levels), set_levels)

  # ── 1. intersection table ──────────────────────────────────────────────────
  inter <- compute_intersections(df, sets)

  inter <- inter |>
    dplyr::filter(count >= min_size) |>
    dplyr::slice_head(n = n_intersections) |>
    dplyr::mutate(intersection_id = dplyr::row_number())

  n_inter <- nrow(inter)
  if (n_inter == 0L) stop("No intersections remain after filtering.")

  # ── 2. set sizes (total membership) ────────────────────────────────────────
  set_sizes <- df |>
    dplyr::select(dplyr::all_of(sets)) |>
    dplyr::summarise(
      dplyr::across(dplyr::everything(), ~ sum(as.logical(.x), na.rm = TRUE))
    ) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to  = "set",
      values_to = "size"
    ) |>
    dplyr::mutate(set = factor(set, levels = set_levels))

  # ── 3. long matrix table for dots & lines ──────────────────────────────────
  matrix_long <- inter |>
    dplyr::select(intersection_id, count, dplyr::all_of(sets)) |>
    tidyr::pivot_longer(
      dplyr::all_of(sets),
      names_to  = "set",
      values_to = "member"
    ) |>
    dplyr::mutate(set = factor(set, levels = set_levels))

  # Segment endpoints derived from the named position lookup, not from
  # factor integer codes, so they are correct regardless of `sets` order.
  segs <- matrix_long |>
    dplyr::filter(member) |>
    dplyr::mutate(y_pos = set_pos[as.character(set)]) |>
    dplyr::group_by(intersection_id) |>
    dplyr::filter(dplyr::n() > 1L) |>
    dplyr::summarise(
      ymin = min(y_pos),
      ymax = max(y_pos),
      .groups = "drop"
    )

  # ── 4. shared x-axis limits ────────────────────────────────────────────────
  x_limits <- c(0.5, n_inter + 0.5)

  # ── 5. TOP panel: intersection size bar chart ──────────────────────────────
  p_top <- ggplot2::ggplot(inter, ggplot2::aes(x = intersection_id, y = count)) +
    ggplot2::geom_col(fill = bar_colour, width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = count),
      vjust = -0.4,
      size  = text_size * 0.25
    ) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(y = "Intersection size", x = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_classic(base_size = text_size) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x  = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(5, 5, 4, 5)
    )

  # ── 6. MATRIX panel: membership dot matrix ─────────────────────────────────
  # Alternating row background rectangles (drawn first, behind dots/lines).
  row_bands <- data.frame(
    y_int = seq_along(set_levels),
    fill  = ifelse(seq_along(set_levels) %% 2 == 0, "#f5f5f5", "white")
  )

  p_matrix <- ggplot2::ggplot(
    matrix_long,
    ggplot2::aes(x = intersection_id, y = set)
  ) +
    # alternating row stripes
    ggplot2::geom_rect(
      data = row_bands,
      ggplot2::aes(
        ymin = y_int - 0.5, ymax = y_int + 0.5,
        xmin = x_limits[1], xmax = x_limits[2],
        fill = I(fill)
      ),
      inherit.aes = FALSE
    ) +
    # absent-set dots
    ggplot2::geom_point(
      data   = dplyr::filter(matrix_long, !member),
      colour = empty_colour,
      size   = 3
    ) +
    # connecting vertical segments (y is numeric position, not factor)
    ggplot2::geom_segment(
      data = segs,
      ggplot2::aes(
        x    = intersection_id, xend = intersection_id,
        y    = ymin,            yend = ymax
      ),
      colour      = dot_colour,
      linewidth   = 1.2,
      inherit.aes = FALSE
    ) +
    # present-set dots
    ggplot2::geom_point(
      data   = dplyr::filter(matrix_long, member),
      colour = dot_colour,
      size   = 3
    ) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(limits = set_levels) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_classic(base_size = text_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      axis.line   = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 5, 5, 5)
    )

  # ── 7. LEFT panel: set size bar chart ──────────────────────────────────────
  p_sets <- ggplot2::ggplot(
    set_sizes,
    ggplot2::aes(x = size, y = set)
  ) +
    ggplot2::geom_col(fill = set_bar_colour, width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = size),
      hjust = -0.2,
      size  = text_size * 0.25
    ) +
    ggplot2::scale_x_reverse(
      expand = ggplot2::expansion(mult = c(0.15, 0))
    ) +
    ggplot2::scale_y_discrete(limits = set_levels) +
    ggplot2::labs(x = "Set size", y = NULL) +
    ggplot2::theme_classic(base_size = text_size) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.line.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(0, 0, 5, 5)
    )

  # ── 8. blank spacer (top-left corner) ──────────────────────────────────────
  p_empty <- ggplot2::ggplot() + ggplot2::theme_void()

  # ── 9. assemble with patchwork ─────────────────────────────────────────────
  # Layout:
  #   [ empty    | top bar chart ]
  #   [ set bars | dot matrix    ]
  (p_empty | p_top) /
    (p_sets | p_matrix) +
    patchwork::plot_layout(
      widths  = c(1, 3),   # set-size bars narrower than main panels
      heights = c(2, 1.2)  # bar chart taller than matrix
    )
}


#' Compute intersection counts from a data frame of logical columns
#'
#' Groups rows by every unique combination of the specified logical columns and
#' returns a summary table with one row per observed intersection, sorted by
#' decreasing count.
#'
#' @param df   A data frame containing (at least) the columns named in
#'   \code{sets}.
#' @param sets Character vector of column names to treat as set-membership
#'   indicators.  Each column must be coercible to \code{logical}.
#'
#' @return A \code{\link[tibble]{tibble}} with one column per set
#'   (\code{TRUE}/\code{FALSE}), a \code{count} column (number of rows in that
#'   intersection), and an \code{intersection_id} integer giving the rank order
#'   (1 = largest intersection).
#' @keywords internal
compute_intersections <- function(df, sets) {
  df |>
    dplyr::select(dplyr::all_of(sets)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.logical)) |>
    dplyr::group_by(dplyr::across(dplyr::everything())) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(count)) |>
    dplyr::mutate(intersection_id = dplyr::row_number())
}

