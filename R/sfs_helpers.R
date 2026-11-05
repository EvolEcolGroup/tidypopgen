# Internal helpers shared by pop_sfs(), pairwise_pop_sfs(), nwise_pop_sfs()
# and pop_sfs_preview().
#
# Missing genotypes are handled as in easySFS (which uses dadi's
# Spectrum.from_data_dict()): each population is projected down to n_proj
# allele copies with the hypergeometric distribution, and loci with fewer
# than n_proj called allele copies in any population of a spectrum are
# dropped from that spectrum. The computation has two stages:
# 1. sfs_config_counts_cpp() counts, in a single pass over the genotypes, the
#    loci with each configuration (number of called and derived allele
#    copies in each population) of each spectrum, where the derived allele is
#    the alternate allele unless the loci are polarised with an ancestral
#    allele (see .sfs_set_ancestral());
# 2. sfs_project_cpp() projects the configurations of each spectrum (in
#    parallel), which are then folded (if needed) in R by .sfs_fold().

#' The genotypes and groups used to compute site frequency spectra
#'
#' @param .x a gen_tibble (grouped or not), or a vctrs_bigSNP
#' @returns a list with `genotypes` (the vctrs_bigSNP), `groups` (the
#'   group information of `.group_info()`, or NULL if not grouped), `ids`
#'   (the 0-based group of each individual), `n_groups`, `labels` (a label for
#'   each group) and `max_n` (the number of allele copies of each group, i.e.
#'   the sum of the ploidy of its individuals), `loci` (the loci table) and
#'   `polarise` (empty, see `.sfs_set_ancestral()`). Empty groups are
#'   dropped, with a message.
#' @keywords internal
#' @noRd
.sfs_data <- function(.x) {
  if (inherits(.x, "grouped_df")) {
    groups <- .group_info(.x)
    # spectra are meant for demographic inference, so (as gt_extract_f2())
    # we drop empty groups (from group_by(.drop = FALSE)), which would be
    # populations without data
    if (any(groups$sizes == 0)) {
      message(
        "empty groups are dropped: ",
        paste(groups$labels[groups$sizes == 0], collapse = ", ")
      )
      .x <- dplyr::group_by(
        .x,
        dplyr::across(dplyr::all_of(dplyr::group_vars(.x))),
        .drop = TRUE
      )
      groups <- .group_info(.x)
    }
    genotypes <- .x$genotypes
    ids <- groups$ids
    labels <- groups$labels
    # the labels name the projections (n_proj), the spectra, the columns of
    # nwise_pop_sfs() and the files of write_sfs(), so they must be distinct
    # (with several grouping variables, their values are pasted with ".", so
    # e.g. ("a.b", "c") and ("a", "b.c") both give "a.b.c")
    dup_labels <- unique(labels[duplicated(labels)])
    if (length(dup_labels) > 0) {
      stop(
        "several groups have the same label (the values of the grouping ",
        "variables pasted with '.'): ", paste(dup_labels, collapse = ", "),
        "; rename the groups"
      )
    }
  } else {
    groups <- NULL
    genotypes <- if (inherits(.x, "vctrs_bigSNP")) .x else .x$genotypes
    ids <- rep(0L, length(genotypes))
    labels <- "all"
  }
  if (length(genotypes) == 0) {
    stop("there are no individuals to compute the site frequency spectrum")
  }
  ploidy <- indiv_ploidy(genotypes)
  # the number of allele copies of each group
  group_of_indiv <- factor(ids, levels = seq_along(labels) - 1L)
  max_n <- as.integer(tapply(ploidy, group_of_indiv, sum))
  names(max_n) <- labels
  list(
    genotypes = genotypes,
    groups = groups,
    ids = ids,
    n_groups = length(labels),
    labels = labels,
    ploidy = ploidy,
    max_n = max_n,
    loci = attr(genotypes, "loci"),
    polarise = integer(0)
  )
}

#' Polarise the loci with their ancestral allele
#'
#' @param data the output of `.sfs_data()`
#' @param ancestral NULL (the reference allele is ancestral), the name of a
#'   column of the loci table with the ancestral allele of each locus, or a
#'   character vector with the ancestral allele of each locus (in the order
#'   of the loci table)
#' @param folded whether the spectra are folded: polarisation does not change
#'   a folded spectrum, so `ancestral` is then ignored (with a warning, as in
#'   dadi, where all loci are used for folded spectra)
#' @returns `data`, with `polarise` set to one value per locus: 0 if the
#'   reference allele is ancestral, 1 if the alternate allele is ancestral,
#'   and -1 if the ancestral allele is missing or is neither allele (the locus
#'   is then dropped, as dadi drops loci without a usable outgroup allele)
#' @keywords internal
#' @noRd
.sfs_set_ancestral <- function(data, ancestral, folded = FALSE) {
  if (is.null(ancestral)) {
    return(data)
  }
  if (folded) {
    warning(
      "`ancestral` is ignored for folded spectra; use `folded = FALSE` to ",
      "polarise the spectra with the ancestral allele"
    )
    return(data)
  }
  loci <- data$loci
  is_column <- is.character(ancestral) && length(ancestral) == 1 &&
    ancestral %in% names(loci)
  if (is_column) {
    ancestral <- loci[[ancestral]]
  }
  if (is.factor(ancestral)) {
    ancestral <- as.character(ancestral)
  }
  if (!is.character(ancestral) && !all(is.na(ancestral))) {
    stop(
      "`ancestral` should be the name of a column of the loci table, or a ",
      "character vector with the ancestral allele of each locus"
    )
  }
  if (length(ancestral) != nrow(loci)) {
    stop(
      "`ancestral` should be the name of a column of the loci table, or have ",
      "one value per locus (", nrow(loci), ")"
    )
  }
  # alleles are compared ignoring case (e.g. lower case ancestral alleles of
  # low confidence)
  ancestral <- toupper(ancestral)
  is_ref <- !is.na(ancestral) & ancestral == toupper(loci$allele_ref)
  is_alt <- !is.na(ancestral) & ancestral == toupper(loci$allele_alt)
  polarise <- rep(-1L, length(ancestral))
  polarise[is_ref %in% TRUE] <- 0L
  polarise[is_alt %in% TRUE] <- 1L
  data$polarise <- polarise
  data
}

#' Count the loci with each configuration of a set of spectra
#'
#' @param data the output of `.sfs_data()`
#' @param spectra a list of integer vectors, each with the (1-based) groups of
#'   a spectrum
#' @param min_n the minimum number of called allele copies, for each group,
#'   for a locus to be counted (loci below it are counted as dropped)
#' @param block_ids optional, the (1-based) bootstrap block of each locus (as
#'   returned by `.jackknife_block_ids()`), to count the configurations in
#'   each block
#' @returns a list with, for each spectrum, `n` and `k` (matrices of the
#'   number of called and alternate allele copies, one row per configuration
#'   and one column per population), `block` (the 0-based block of each
#'   configuration), `count` (the number of loci with each configuration),
#'   and `n_dropped` and `n_unpolarised` (the number of loci dropped in each
#'   block)
#' @keywords internal
#' @noRd
.sfs_config_counts <- function(data, spectra, min_n, n_cores, block_size,
                               block_ids = NULL) {
  genotypes <- data$genotypes
  if (is.null(block_ids)) {
    block_ids <- integer(0)
    n_blocks <- 1L
  } else {
    n_blocks <- max(block_ids)
    block_ids <- as.integer(block_ids) - 1L
  }
  sfs_config_counts_cpp(
    BM = attr(genotypes, "fbm"),
    rowInd = vctrs::vec_data(genotypes),
    colInd = attr(genotypes, "loci")$big_index,
    groupIds = data$ids,
    ngroups = data$n_groups,
    ploidy = data$ploidy,
    is_pseudohaploid = is_pseudohaploid(genotypes),
    spectra = lapply(spectra, function(x) as.integer(x) - 1L),
    min_n = as.integer(min_n),
    polarise = data$polarise,
    block_ids = block_ids,
    n_blocks = as.integer(n_blocks),
    ncores = n_cores,
    # the kernel keeps the counts of a block of loci for all groups, so we
    # cap the size of the blocks
    block_size = as.integer(min(block_size, 4096))
  )
}

#' Check the number of allele copies to project each group to
#'
#' @param n_proj a single value (used for all groups), an unnamed vector with
#'   one value per group (in the order of the groups), or a vector named with
#'   the group labels
#' @param max_n the number of allele copies of each group (named with the
#'   group labels)
#' @returns an integer vector with one value per group, named with the group
#'   labels
#' @keywords internal
#' @noRd
.sfs_check_n_proj <- function(n_proj, max_n) {
  labels <- names(max_n)
  is_whole <- is.numeric(n_proj) && length(n_proj) > 0 && !anyNA(n_proj) &&
    all(n_proj == round(n_proj))
  if (!is_whole) {
    stop("`n_proj` should be a vector of whole numbers")
  }
  if (!is.null(names(n_proj))) {
    unknown <- setdiff(names(n_proj), labels)
    if (length(unknown) > 0) {
      stop(
        "`n_proj` has names that are not groups: ",
        paste(unknown, collapse = ", ")
      )
    }
    missing_groups <- setdiff(labels, names(n_proj))
    if (length(missing_groups) > 0) {
      stop(
        "`n_proj` has no value for: ",
        paste(missing_groups, collapse = ", ")
      )
    }
    n_proj <- n_proj[labels]
  } else if (length(n_proj) == 1) {
    n_proj <- rep(n_proj, length(labels))
  } else if (length(n_proj) != length(labels)) {
    stop("`n_proj` should have a single value, or one value per group")
  }
  if (any(n_proj < 1)) {
    stop("`n_proj` should be at least 1")
  }
  too_large <- n_proj > max_n
  if (any(too_large)) {
    stop(
      "`n_proj` is larger than the number of allele copies of: ",
      paste0(labels[too_large], " (", max_n[too_large], ")", collapse = ", ")
    )
  }
  n_proj <- as.integer(n_proj)
  names(n_proj) <- labels
  n_proj
}

#' Preview the projections of each group (easySFS --preview)
#'
#' @param configs the configurations of the 1D spectrum of each group (as
#'   returned by `.sfs_config_counts()`, with `min_n = 1`)
#' @param max_n the number of allele copies of each group
#' @returns a list with, for each group, a data.frame with columns `n_proj`
#'   (from 2 to the number of allele copies of the group), `n_seg_sites` and
#'   `n_loci`, and a logical column `best` (see `.sfs_best_proj()`)
#' @keywords internal
#' @noRd
.sfs_preview_from_configs <- function(configs, max_n) {
  lapply(seq_along(configs), function(i_grp) {
    conf <- configs[[i_grp]]
    preview <- sfs_preview_cpp(
      n = conf$n[, 1],
      k = conf$k[, 1],
      count = conf$count,
      n_max = max_n[i_grp]
    )
    # projections from 2 copies (as easySFS)
    p <- seq_len(max_n[i_grp])[-1]
    res <- data.frame(
      n_proj = p,
      n_seg_sites = preview$n_seg_sites[p + 1],
      n_loci = preview$n_loci[p + 1]
    )
    res$best <- res$n_proj == .sfs_best_proj(res, max_n[i_grp])
    res
  })
}

#' Choose the number of allele copies to project a group to
#'
#' Following the rationale of easySFS and of the dadi manual, we choose the
#' projection that maximises the number of segregating sites in the projected
#' spectrum: projecting to more allele copies uncovers more rare variants,
#' but drops more loci with missing data. Ties are broken by choosing the
#' largest projection (the largest sample size).
#'
#' @param preview a data.frame with columns `n_proj` and `n_seg_sites`
#' @param max_n the number of allele copies of the group (used if there are
#'   fewer than 2, i.e. when preview has no rows)
#' @returns a single integer
#' @keywords internal
#' @noRd
.sfs_best_proj <- function(preview, max_n) {
  if (nrow(preview) == 0) {
    return(as.integer(max_n))
  }
  is_max <- preview$n_seg_sites == max(preview$n_seg_sites)
  as.integer(max(preview$n_proj[is_max]))
}

#' The projection of each group
#'
#' If `n_proj` is NULL, it is chosen for each group with `.sfs_best_proj()`
#' (which requires a pass over the genotypes, unless the configurations of
#' the 1D spectra are given), and a message lists the chosen values.
#'
#' @param n_proj NULL, or the projections given by the user
#' @param data the output of `.sfs_data()`
#' @param configs_1d optional, the configurations of the 1D spectrum of each
#'   group (computed with `min_n = 1`)
#' @returns a named integer vector (see `.sfs_check_n_proj()`)
#' @keywords internal
#' @noRd
.sfs_n_proj <- function(n_proj, data, n_cores, block_size,
                        configs_1d = NULL) {
  if (!is.null(n_proj)) {
    return(.sfs_check_n_proj(n_proj, data$max_n))
  }
  if (is.null(configs_1d)) {
    configs_1d <- .sfs_config_counts(data,
      spectra = as.list(seq_len(data$n_groups)),
      min_n = rep(1L, data$n_groups),
      n_cores = n_cores,
      block_size = block_size
    )
  }
  previews <- .sfs_preview_from_configs(configs_1d, data$max_n)
  n_proj <- vapply(seq_along(previews), function(i_grp) {
    .sfs_best_proj(previews[[i_grp]], data$max_n[i_grp])
  }, integer(1))
  names(n_proj) <- data$labels
  message(
    "`n_proj` not given; each group is projected to the number of allele ",
    "copies that maximises its number of segregating sites (see ",
    "pop_sfs_preview()): ",
    paste(names(n_proj), n_proj, sep = " = ", collapse = ", ")
  )
  n_proj
}

#' Project the configurations of a set of spectra, for each replicate
#'
#' @param configs the configurations of each spectrum (as returned by
#'   `.sfs_config_counts()`)
#' @param spectra a list of integer vectors, each with the (1-based) groups of
#'   a spectrum
#' @param n_proj the projection of each group (named with the group labels)
#' @param folded whether to fold the spectra
#' @param n_cores the number of cores used to project the spectra
#' @param weights a matrix with the weight of each block (columns) in each
#'   bootstrap replicate (rows); by default, a single replicate with all loci
#'   (a single block of weight 1)
#' @returns a list with, for each spectrum, `values` (a matrix of the cells of
#'   the spectrum, in column-major order, for each replicate (columns)),
#'   `mask`, `dims`, `pop_ids`, `folded`, and the number of loci used
#'   (`n_loci`), dropped for missing data (`n_loci_dropped`) and without
#'   ancestral allele (`n_loci_unpolarised`) in each replicate
#' @keywords internal
#' @noRd
.sfs_project <- function(configs, spectra, n_proj, folded, n_cores,
                         weights = matrix(1, 1, 1)) {
  storage.mode(weights) <- "double"
  projected <- sfs_project_cpp(
    configs = configs,
    spectra = lapply(spectra, function(x) as.integer(x) - 1L),
    n_proj = unname(n_proj),
    weights = weights,
    ncores = n_cores
  )
  lapply(seq_along(spectra), function(i) {
    this_proj <- n_proj[spectra[[i]]]
    dims <- unname(this_proj) + 1L
    values <- projected[[i]]$sfs
    # mask the corners (monomorphic in the projected sample), as dadi does
    mask <- array(FALSE, dim = dims)
    mask[1] <- TRUE
    mask[length(mask)] <- TRUE
    if (folded) {
      folded_values <- .sfs_fold_values(values, mask)
      values <- folded_values$values
      mask <- folded_values$mask
    }
    # the loci of each replicate are the loci of its blocks
    block_sums <- function(x) as.vector(weights %*% x)
    list(
      values = values,
      mask = mask,
      dims = dims,
      pop_ids = names(this_proj),
      folded = folded,
      n_loci = block_sums(projected[[i]]$n_used),
      n_loci_dropped = block_sums(
        projected[[i]]$n_dropped + configs[[i]]$n_dropped
      ),
      n_loci_unpolarised = block_sums(configs[[i]]$n_unpolarised)
    )
  })
}

#' The `sfs` object of a replicate of a projected spectrum
#'
#' @param projected an element of the list returned by `.sfs_project()`
#' @param i_rep the replicate
#' @returns an `sfs` object
#' @keywords internal
#' @noRd
.sfs_replicate <- function(projected, i_rep = 1) {
  new_sfs(
    array(projected$values[, i_rep], dim = projected$dims),
    mask = projected$mask,
    pop_ids = projected$pop_ids,
    folded = projected$folded,
    n_loci = projected$n_loci[i_rep],
    n_loci_dropped = projected$n_loci_dropped[i_rep],
    n_loci_unpolarised = projected$n_loci_unpolarised[i_rep]
  )
}

#' Build spectra from their configurations
#'
#' @inheritParams .sfs_project
#' @returns a list of `sfs` objects (computed with all loci)
#' @keywords internal
#' @noRd
.sfs_from_configs <- function(configs, spectra, n_proj, folded, n_cores) {
  lapply(
    .sfs_project(configs, spectra, n_proj, folded, n_cores),
    .sfs_replicate
  )
}

#' Fold a spectrum
#'
#' The algorithm of dadi's `Spectrum.fold()`: a joint spectrum is folded on
#' the total number of (alternate) allele copies across populations. Cells
#' with more than half of the total number of copies are added to their
#' mirror cell (the cell obtained by reversing the array along all axes),
#' and set to 0; cells with exactly half the copies are averaged with their
#' mirror cell. The mask covers the cells masked before folding, their
#' mirror cells, and the folded out cells.
#'
#' @param fs a numeric array (one dimension per population)
#' @param mask a logical array of the same dimensions
#' @returns a list with the folded `fs` and its `mask`
#' @keywords internal
#' @noRd
.sfs_fold <- function(fs, mask) {
  folded <- .sfs_fold_values(matrix(as.vector(fs), ncol = 1), mask)
  list(
    fs = array(folded$values[, 1], dim = dim(mask)),
    mask = folded$mask
  )
}

#' Fold the spectra of several replicates
#'
#' As `.sfs_fold()`, for a matrix with the cells of a spectrum (in
#' column-major order) for each replicate (columns).
#'
#' @param values a numeric matrix, with one row per cell and one column per
#'   replicate
#' @param mask a logical array with the dimensions of the spectrum
#' @returns a list with the folded `values` and the `mask`
#' @keywords internal
#' @noRd
.sfs_fold_values <- function(values, mask) {
  dims <- dim(mask)
  total_samples <- sum(dims - 1)
  # the total number of copies of each cell
  total_per_entry <- array(0L, dim = dims)
  for (i in seq_along(dims)) {
    total_per_entry <- total_per_entry + slice.index(mask, i) - 1L
  }
  # reversing an array along all its axes reverses its elements
  rev_idx <- rev(seq_along(mask))
  where_folded_out <- as.vector(total_per_entry > floor(total_samples / 2))
  folded <- values + (values * where_folded_out)[rev_idx, , drop = FALSE]
  folded[where_folded_out, ] <- 0
  where_ambiguous <- as.vector(total_per_entry == total_samples / 2)
  ambiguous <- values * where_ambiguous
  folded <- folded - 0.5 * ambiguous + 0.5 * ambiguous[rev_idx, , drop = FALSE]
  list(
    values = folded,
    mask = array(mask | rev(mask) | where_folded_out, dim = dims)
  )
}

#' Create an `sfs` object
#'
#' @param fs a numeric array, with one dimension per population (of extent
#'   the projection plus 1)
#' @param mask a logical array of the same dimensions, TRUE for the masked
#'   cells (the corners and, for folded spectra, the folded out cells)
#' @param pop_ids the population labels
#' @param folded whether the spectrum is folded
#' @param n_loci the number of loci used (projected)
#' @param n_loci_dropped the number of loci dropped (too few called copies in
#'   at least one population)
#' @param n_loci_unpolarised the number of loci dropped because their
#'   ancestral allele is unknown (see `.sfs_set_ancestral()`)
#' @returns an `sfs` object: the array, with dimnames giving the number of
#'   allele copies in each population, and the other arguments as attributes
#'   (plus `n_proj`)
#' @keywords internal
#' @noRd
new_sfs <- function(fs, mask, pop_ids, folded, n_loci, n_loci_dropped,
                    n_loci_unpolarised = 0) {
  dims <- dim(fs)
  dimnames(fs) <- stats::setNames(
    lapply(dims, function(d) as.character(seq_len(d) - 1L)),
    pop_ids
  )
  structure(
    fs,
    mask = unname(mask),
    pop_ids = pop_ids,
    n_proj = stats::setNames(as.integer(dims - 1L), pop_ids),
    folded = folded,
    n_loci = n_loci,
    n_loci_dropped = n_loci_dropped,
    n_loci_unpolarised = n_loci_unpolarised,
    class = "sfs"
  )
}

#' The order of the cells of a spectrum in a tibble
#'
#' @param dims the dimensions of the spectrum
#' @param fastest the population whose number of allele copies changes
#'   fastest from one row to the next: "last" (the order of the populations
#'   is respected: rows are sorted by the first population, then the second,
#'   ...), or "first"
#' @param count_names the names of the columns with the number of allele
#'   copies of each population
#' @returns a list with `counts` (a tibble with the number of copies of each
#'   population in each cell, one row per cell) and `order` (the index, in
#'   the column-major array, of the cell of each row)
#' @keywords internal
#' @noRd
.sfs_cell_index <- function(dims, fastest = c("last", "first"), count_names) {
  fastest <- match.arg(fastest)
  copies <- lapply(dims, function(d) seq_len(d) - 1L)
  if (fastest == "first") {
    # R arrays are column-major: the first dimension changes fastest
    counts <- expand.grid(copies, KEEP.OUT.ATTRS = FALSE)
    cell_order <- seq_len(prod(dims))
  } else {
    perm <- rev(seq_along(dims))
    counts <- expand.grid(copies[perm], KEEP.OUT.ATTRS = FALSE)[perm]
    cell_order <- as.vector(aperm(array(seq_len(prod(dims)), dim = dims), perm))
  }
  names(counts) <- count_names
  list(counts = tibble::as_tibble(counts), order = cell_order)
}

#' The cells of a spectrum as a tibble
#'
#' @param x an `sfs` object
#' @inheritParams .sfs_cell_index
#' @returns a tibble with one row per cell, one column of allele copies per
#'   population, `n_sites` (the projected number of sites) and `masked`
#' @keywords internal
#' @noRd
.sfs_cells <- function(x, fastest = c("last", "first"), count_names) {
  index <- .sfs_cell_index(dim(x), fastest, count_names)
  index$counts %>%
    dplyr::mutate(
      n_sites = as.vector(unclass(x))[index$order],
      masked = as.vector(attr(x, "mask"))[index$order]
    )
}

#' The cells of the bootstrap replicates of a spectrum as a tibble
#'
#' @param projected an element of the list returned by `.sfs_project()`
#' @inheritParams .sfs_cell_index
#' @returns a tibble with one row per replicate and cell, with columns
#'   `replicate`, one column of allele copies per population, `n_sites` and
#'   `masked`
#' @keywords internal
#' @noRd
.sfs_replicate_cells <- function(projected, fastest = c("last", "first"),
                                 count_names) {
  index <- .sfs_cell_index(projected$dims, fastest, count_names)
  n_rep <- ncol(projected$values)
  n_cells <- length(index$order)
  rows <- rep(seq_len(n_cells), times = n_rep)
  dplyr::bind_cols(
    tibble::tibble(replicate = rep(seq_len(n_rep), each = n_cells)),
    index$counts[rows, , drop = FALSE],
    tibble::tibble(
      n_sites = as.vector(projected$values[index$order, , drop = FALSE]),
      masked = rep(as.vector(projected$mask)[index$order], times = n_rep)
    )
  )
}

#' Tidy a site frequency spectrum
#'
#' Turn an `sfs` object (as returned by [pop_sfs()], [pairwise_pop_sfs()] and
#' [nwise_pop_sfs()] with `type = "array"`) into a tibble with one row per
#' cell of the spectrum.
#'
#' @param x an `sfs` object
#' @param fastest which population's number of allele copies changes fastest
#'   from one row to the next: "last" (the default, which keeps the order of
#'   the populations: rows are sorted by the number of copies in the first
#'   population, then in the second, and so on, as in the files written by
#'   dadi and easySFS) or "first".
#' @param ... currently unused.
#' @returns a tibble with one column per population (named with the
#'   population labels) giving the number of allele copies (alternate copies
#'   for an unfolded spectrum, copies of the minor allele across all the
#'   populations of the spectrum for a folded one), `n_sites`, the number
#'   of sites in that cell (a fractional number, as each locus is spread over
#'   the cells by the projection), and `masked`, TRUE for the cells that are
#'   masked (the corners, i.e. monomorphic sites, and for a folded spectrum
#'   the cells that were folded out). The projection of each population is
#'   stored in the attribute `n_proj`, and whether the spectrum is folded in
#'   the attribute `folded`.
#' @export
#' @examples
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' sfs <- example_gt %>% pop_sfs(n_proj = 4, type = "array")
#' tidy(sfs[[1]])
tidy.sfs <- function(x, fastest = c("last", "first"), ...) {
  rlang::check_dots_empty()
  .sfs_tidy_attributes(
    .sfs_cells(x, fastest = fastest, count_names = attr(x, "pop_ids")),
    n_proj = attr(x, "n_proj"),
    folded = attr(x, "folded")
  )
}

#' @export
print.sfs <- function(x, ...) {
  n_proj <- attr(x, "n_proj")
  cat(paste0(
    "# ", if (attr(x, "folded")) "A folded" else "An unfolded",
    " site frequency spectrum of ", length(n_proj),
    if (length(n_proj) == 1) " population" else " populations", "\n",
    "# Allele copies: ",
    paste(names(n_proj), n_proj, sep = " = ", collapse = ", "), "\n",
    "# Loci used: ", format(attr(x, "n_loci")),
    " (dropped for missing data: ", format(attr(x, "n_loci_dropped")),
    if (isTRUE(attr(x, "n_loci_unpolarised") > 0)) {
      paste0(
        "; without ancestral allele: ",
        format(attr(x, "n_loci_unpolarised"))
      )
    },
    ")\n"
  ))
  y <- unclass(x)
  attributes(y) <- list(dim = dim(x), dimnames = dimnames(x))
  print(y, ...)
  invisible(x)
}
