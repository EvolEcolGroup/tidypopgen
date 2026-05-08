# UpSet Plot using ggplot2 + patchwork
# Inputs: a data frame with logical columns representing set membership
# No specialised packages (no UpSetR, no ComplexUpset)

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# ── helpers ──────────────────────────────────────────────────────────────────

#' Compute all intersection counts from a data frame of logical columns.
#'
#' @param df   Data frame containing (at least) the columns named in `sets`.
#' @param sets Character vector of logical column names (the sets).
#' @return A tibble with one row per observed intersection containing:
#'   - one column per set (TRUE/FALSE membership),
#'   - `count`      : number of rows in that intersection,
#'   - `intersection_id` : pasted set-membership string (for ordering).
compute_intersections <- function(df, sets) {
  df |>
    select(all_of(sets)) |>
    mutate(across(everything(), as.logical)) |>
    group_by(across(everything())) |>
    summarise(count = n(), .groups = "drop") |>
    arrange(desc(count)) |>
    mutate(intersection_id = row_number())
}

#' Build the three ggplot panels and combine with patchwork.
#'
#' @param df            Data frame with logical set-membership columns.
#' @param sets          Character vector of column names to use as sets.
#'                      If NULL all logical columns are used.
#' @param min_size      Drop intersections with fewer than this many members.
#' @param n_intersections Keep only the top-n intersections (by count).
#' @param bar_colour    Fill colour for the intersection bar chart.
#' @param dot_colour    Colour for filled dots / connecting lines in the matrix.
#' @param empty_colour  Colour for unfilled (absent) dots in the matrix.
#' @param set_bar_colour Fill colour for the set-size bar chart.
#' @param text_size     Base text size passed to ggplot themes.
#' @return A patchwork object that can be printed or saved with ggsave().
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

  # ── 0. resolve set columns ────────────────────────────────────────────────
  if (is.null(sets)) {
    sets <- names(df)[sapply(df, is.logical)]
    if (length(sets) == 0L)
      stop("No logical columns found. Supply `sets` explicitly.")
  }

  # ── 1. intersection table ─────────────────────────────────────────────────
  inter <- compute_intersections(df, sets)

  # filter & truncate
  inter <- inter |>
    filter(count >= min_size) |>
    slice_head(n = n_intersections) |>
    mutate(intersection_id = row_number())   # re-index after filtering

  n_inter <- nrow(inter)
  if (n_inter == 0L) stop("No intersections remain after filtering.")

  # ── 2. set sizes (total membership) ───────────────────────────────────────
  set_sizes <- df |>
    select(all_of(sets)) |>
    summarise(across(everything(), ~ sum(as.logical(.x), na.rm = TRUE))) |>
    pivot_longer(everything(), names_to = "set", values_to = "size") |>
    mutate(set = factor(set, levels = rev(sets)))   # rev so top set is first

  # ── 3. long matrix table for dots & lines ────────────────────────────────
  matrix_long <- inter |>
    select(intersection_id, count, all_of(sets)) |>
    pivot_longer(all_of(sets), names_to = "set", values_to = "member") |>
    mutate(set = factor(set, levels = rev(sets)))

  # segment endpoints: connect dots that are TRUE within each intersection
  segs <- matrix_long |>
    filter(member) |>
    group_by(intersection_id) |>
    filter(n() > 1) |>          # only draw lines when ≥2 sets active
    summarise(
      ymin = min(as.integer(set)),
      ymax = max(as.integer(set)),
      .groups = "drop"
    )

  # ── 4. shared x-axis limits ───────────────────────────────────────────────
  x_limits <- c(0.5, n_inter + 0.5)

  # ── 5. TOP panel: intersection size bar chart ─────────────────────────────
  p_top <- ggplot(inter, aes(x = intersection_id, y = count)) +
    geom_col(fill = bar_colour, width = 0.6) +
    geom_text(aes(label = count), vjust = -0.4, size = text_size * 0.25) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(y = "Intersection size", x = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = text_size) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank(),
      plot.margin  = margin(5, 5, 4, 5)
    )

  # -- 6. MIDDLE panel: membership matrix ----------------------------------
  # alternating row background rectangles
  set_levels <- levels(matrix_long$set)
  row_bands <- data.frame(
    y_int = seq_along(set_levels),
    fill  = ifelse(seq_along(set_levels) %% 2 == 0, "#f5f5f5", "white")
  )

  p_matrix <- ggplot(matrix_long, aes(x = intersection_id, y = set)) +
    # alternating row stripes (drawn first, behind everything)
    geom_rect(
      data = row_bands,
      aes(ymin = y_int - 0.5, ymax = y_int + 0.5,
          xmin = x_limits[1],  xmax = x_limits[2],
          fill = I(fill)),
      inherit.aes = FALSE
    ) +
    # background grid dots (absent)
    geom_point(
      data = filter(matrix_long, !member),
      colour = empty_colour, size = 3
    ) +
    # connecting vertical segments
    geom_segment(
      data = segs,
      aes(x = intersection_id, xend = intersection_id,
          y = ymin, yend = ymax),
      colour = dot_colour, linewidth = 1,
      inherit.aes = FALSE
    ) +
    # filled dots (present)
    geom_point(
      data = filter(matrix_long, member),
      colour = dot_colour, size = 3
    ) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_discrete() +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = text_size) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks   = element_blank(),
      axis.line    = element_blank(),
      plot.margin  = margin(0, 5, 5, 5)
    )


  # ── 7. RIGHT panel: set size bar chart ────────────────────────────────────
  p_sets <- ggplot(set_sizes, aes(x = size, y = set)) +
    geom_col(fill = set_bar_colour, width = 0.6) +
#    geom_text(aes(label = size), hjust = -0.2, size = text_size * 0.25) +
    scale_x_reverse(expand = expansion(mult = c(0.15, 0))) +
    labs(x = "Set size", y = NULL) +
    theme_classic(base_size = text_size) +
    theme(
      axis.text.y  = element_blank(),
      axis.line.y  = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin  = margin(0, 0, 5, 5)
    )

  # ── 8. empty spacer (top-left corner) ─────────────────────────────────────
  p_empty <- ggplot() + theme_void()

  # ── 9. assemble with patchwork ────────────────────────────────────────────
  # Layout:
  #   [ empty  |   top bar chart  ]
  #   [ set bars | matrix         ]
  (p_empty | p_top) /
    (p_sets | p_matrix) +
    plot_layout(
      widths  = c(1, 3),   # set-bars narrower than main panels
      heights = c(2, 1.2)  # bar chart taller than matrix
    )
}


# ── demo ─────────────────────────────────────────────────────────────────────

set.seed(42)
n <- 200

demo_df <- tibble(
  SetA = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.6, 0.4)),
  SetB = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.5, 0.5)),
  SetC = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.4, 0.6)),
  SetD = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.3, 0.7)),
  extra_col = rnorm(n)   # non-logical column — ignored automatically
)

p <- upset_plot(demo_df, text_size = 10)
print(p)
