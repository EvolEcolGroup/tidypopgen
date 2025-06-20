#' Autoplots for `gt_dapc` objects
#'
#' For `gt_dapc`, the following types of plots are available:
#' - `screeplot`: a plot of the eigenvalues of the discriminant axes
#' - `scores` a scatterplot of the scores of each individual on two discriminant
#' axes (defined by `ld`)
#' - `loadings` a plot of loadings of all loci for a discriminant axis
#'    (chosen with `ld`)
#' - `components` a bar plot showing the probability of assignment to
#'    each cluster
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are not
#' customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' @param object an object of class `gt_dapc`
#' @param type the type of plot (one of "screeplot", "scores", "loadings", and
#'  "components")
#' @param ld the principal components to be plotted: for scores, a pair of
#'   values e.g. c(1,2); for `loadings` either one or more values.
#' @param group a vector of group memberships to order the individuals in
#'   "components" plot. If NULL, the clusters used for the DAPC will be used.
#' @param n_col for `loadings` plots, if multiple LD axis are plotted, how many
#'   columns should be used.
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @rdname autoplot_gt_dapc
#' @export
#' @examples
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' lobsters <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Remove monomorphic loci and impute
#' lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
#' lobsters <- gt_impute_simple(lobsters, method = "mode")
#'
#' # Create PCA and run DAPC
#' pca <- gt_pca_partialSVD(lobsters)
#' populations <- as.factor(lobsters$population)
#' dapc_res <- gt_dapc(pca, n_pca = 6, n_da = 2, pop = populations)
#'
#' # Screeplot
#' autoplot(dapc_res, type = "screeplot")
#'
#' # Scores plot
#' autoplot(dapc_res, type = "scores", ld = c(1, 2))
#'
#' # Loadings plot
#' autoplot(dapc_res, type = "loadings", ld = 1)
#'
#' # Components plot
#' autoplot(dapc_res, type = "components", group = populations)
#'
autoplot.gt_dapc <- function(
    object,
    type = c(
      "screeplot",
      "scores",
      "loadings",
      "components"
    ),
    ld = NULL,
    group = NULL,
    n_col = 1,
    ...) {
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (type == "screeplot") {
    tidy(object, matrix = "eigenvalues") %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$LD, y = .data$eigenvalue)) +
      ggplot2::geom_point() +
      ggplot2::geom_line()
  } else if (type == "scores") {
    if (is.null(ld)) {
      ld <- c(1, 2)
    }
    if (length(ld) != 2) {
      stop("for 'scores' plots, 'ld' should be a pair of values, e.g. c(1,2)")
    }
    tibble(cluster = object$grp) %>%
      mutate(
        LDa = object$ind.coord[, ld[1]],
        LDb = object$ind.coord[, ld[2]]
      ) %>%
      ggplot2::ggplot(ggplot2::aes(
        x = .data$LDa,
        y = .data$LDb,
        colour = .data$cluster
      )) +
      ggplot2::geom_point() +
      ggplot2::stat_ellipse() +
      ggplot2::labs(x = paste0("LD", ld[1]), y = paste0("LD", ld[2]))
  } else if (type == "loadings") {
    # code modified from bigstatsr for consistency with pca
    # stop("autoplot for gt_dapc does not have a loading option yet") # nolint
    if (is.null(object$var.load)) {
      stop(paste(
        "the dapc object was saved without loadings for loci,",
        "rerun it with 'loadings_by_locus = TRUE'"
      ))
    }
    if (is.null(ld)) {
      ld <- 1
    }
    if (length(ld) > 1) {
      all.p <- lapply(ld, function(i) {
        p <- autoplot(object, type = "loadings", ld = i)
        p$layers[[1]] <- NULL
        p + ggplot2::geom_hex() + ggplot2::scale_fill_viridis_c()
      })

      patchwork::wrap_plots(all.p, ncol = n_col)
    } else {
      p <- ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = bigstatsr::rows_along(object$var.load),
          y = object$var.load[, ld]
        )
      ) +
        ggplot2::geom_point() +
        bigstatsr::theme_bigstatsr(size.rel = 1) +
        ggplot2::labs(
          title = paste0("Loadings of LD", ld),
          x = "Locus",
          y = NULL
        )
      nval <- nrow(object$var.load)
      `if`(
        nval > 12,
        p,
        p + ggplot2::scale_x_discrete(limits = factor(seq_len(nval)))
      )
    }
  } else if (type == "components") {
    if (is.null(group)) {
      group <- object$grp
    }
    q_tbl <- object$posterior %>%
      tibble::as_tibble() %>%
      dplyr::rename_with(~ paste0("Q", .x)) %>%
      # add the pops data for plotting
      dplyr::mutate(
        name = rownames(object$posterior),
        group = group
      ) %>%
      tidyr::pivot_longer(
        cols = starts_with("Q"),
        names_to = "q",
        values_to = "prob"
      )

    ggplot2::ggplot(
      q_tbl,
      ggplot2::aes(
        x = .data$name,
        y = .data$prob,
        fill = .data$q
      )
    ) +
      ggplot2::geom_col(color = "gray", size = 0.1) +
      ggplot2::facet_grid(
        ~group,
        switch = "x",
        scales = "free",
        space = "free"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "", y = "") +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.7)) +
      ggplot2::theme(
        panel.spacing.x = ggplot2::unit(0.01, "lines"),
        axis.text.x = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::guides(fill = "none")
  }
}
