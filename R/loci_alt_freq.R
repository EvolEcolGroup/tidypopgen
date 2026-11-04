#' Estimate allele frequencies at each locus
#'
#' Allele frequencies can be estimates as minimum allele frequencies (MAF) with
#' `loci_maf()` or the frequency of the alternate allele (with
#' `loci_alt_freq()`). The latter are in line with the genotypes matrix (e.g. as
#' extracted by [`show_loci()`]). Most users will be in interested in the MAF,
#' but the raw frequencies might be useful when computing aggregated statistics.
#' Both `loci_maf()` and `loci_alt_freq()` have efficient methods to support
#' grouped `gen_tibble` objects. These can return a tidied tibble, a list, or a
#' matrix.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#'   directly to the function). This defaults to "genotypes" and can only take
#'   that value. There is no need for the user to set it, but it is included to
#'   resolve certain tidyselect operations.
#' @param as_counts boolean defining whether the count of alternate and valid
#'   (i.e. total number) alleles (rather than the frequencies) should be
#'   returned. It defaults to FALSE (i.e. frequencies are returned by default).
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param block_size maximum number of loci read at once.
#' @param type type of object to return, if using grouped method. One of "tidy",
#'   "list", or "matrix". Default is "tidy".
#'   Data can be grouped by more than one variable: groups are then labelled
#'   by the values of the grouping variables pasted with "." (e.g. `pop1.F`),
#'   and the "tidy" output has one column per grouping variable rather than a
#'   single `group` column. The `group` column, or the columns of the
#'   grouping variables, keep the type of the grouping variables (e.g. a
#'   factor keeps its levels), so that groups stay in their order when the
#'   output is sorted or plotted.
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns a vector of frequencies, one per locus, if `as_counts = FALSE`;
#' else a matrix of two columns, the count of alternate alleles and the count
#' valid alleles (i.e. the sum of alternate and reference)
#' @rdname loci_alt_freq
#' @export
#' @examples
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # For alternate allele frequency
#' example_gt %>% loci_alt_freq()
#'
#' # For alternate allele frequency per locus per population
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_alt_freq()
#' # alternatively, return a list of populations with their frequencies
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_alt_freq(type = "list")
#' # or a matrix with populations in columns and loci in rows
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_alt_freq(type = "matrix")
#' # or within reframe (not recommended, as it much less efficient
#' # than using it directly as shown above)
#' library(dplyr)
#' example_gt %>%
#'   group_by(population) %>%
#'   reframe(alt_freq = loci_alt_freq(genotypes))

#' # For MAF
#' example_gt %>% loci_maf()
#'
#' # For minor allele frequency per locus per population
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_maf()
#' # alternatively, return a list of populations with their frequencies
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_maf(type = "list")
#' # or a matrix with populations in columns and loci in rows
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_maf(type = "matrix")
#'
loci_alt_freq <- function(.x,
                          .col = "genotypes",
                          as_counts = FALSE,
                          n_cores, block_size, type, ...) {
  UseMethod("loci_alt_freq", .x)
}

#' @export
#' @rdname loci_alt_freq
loci_alt_freq.tbl_df <- function(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  # multicore is used by openMP within the
  # freq cpp function
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  # the bigapply that splits in blocks is not
  # multithreaded, as we use the multiple threads
  # for openMP
  ...
) {
  stopifnot_gen_tibble(.x)
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_alt_freq only works with the genotypes column")
  }
  loci_alt_freq(.x$genotypes,
    as_counts = as_counts,
    n_cores = n_cores, block_size = block_size
  )
}


#' @export
#' @rdname loci_alt_freq
loci_alt_freq.vctrs_bigSNP <- function(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
) {
  rlang::check_dots_empty()
  # if we have diploid
  if (is_diploid_only(.x) || is_pseudohaploid(.x)) {
    ploidy <- indiv_ploidy(.x)
    loci_alt_freq_dip_pseudo(.x,
      n_cores = n_cores, block_size = block_size,
      ploidy = ploidy,
      is_pseudohaploid = is_pseudohaploid(.x),
      as_counts = as_counts
    )
  } else {
    loci_alt_freq_polyploid(.x,
      as_counts = as_counts,
      n_cores = n_cores, block_size = block_size, ...
    )
  }
}

#' @export
#' @rdname loci_alt_freq
loci_alt_freq.grouped_df <- function(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  type = c("tidy", "list", "matrix"),
  ...
) {
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_alt_freq only works with the genotypes column")
  }

  rlang::check_dots_empty()
  type <- match.arg(type)
  if (is_diploid_only(.x) || is_pseudohaploid(.x)) {
    # number of groups (used to define dimensions of objects)
    n_groups <- dplyr::n_groups(.x)

    freq_mat <- grouped_alt_freq_dip_pseudo(
      .x,
      as_counts = as_counts,
      n_cores = n_cores,
      block_size = block_size
    )

    if (as_counts) {
      # split into two matrices (alt counts and valid allele counts, and return
      # as a list)
      counts_list <- list(
        n_alt = freq_mat[, 1:n_groups],
        n_valid = freq_mat[, (n_groups + 1):(n_groups * 2)]
      )
      # add col names
      names(counts_list) <- c("n_alt", "n_valid")
      return(counts_list)
    }

    groups <- .group_info(.x)
    freq_mat <- format_grouped_output(
      out_mat = freq_mat,
      group_ids = groups$labels,
      loci_names = loci_names(.x),
      type = type,
      group_keys = groups$keys
    )
  } else {
    # the polyploid case
    stop(
      "loci_alt_freq for polyploid is not directly implemented yet,",
      "use group_map(.x, .f = ~ loci_alt_freq(.x)) for polyploid cases"
    )
  }
  return(freq_mat)
}

#' @rdname loci_alt_freq
#' @export
loci_maf <- function(.x,
                     .col = "genotypes",
                     n_cores, block_size, type, ...) {
  UseMethod("loci_maf", .x)
}

#' @export
#' @rdname loci_alt_freq
loci_maf.tbl_df <- function(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
) {
  stopifnot_gen_tibble(.x) # confirm that .col is "genotypes"
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_maf only works with the genotypes column")
  }
  loci_maf(.x$genotypes, n_cores = n_cores, block_size = block_size, ...)
}

#' @export
#' @rdname loci_alt_freq
loci_maf.vctrs_bigSNP <- function(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
) {
  freq <- loci_alt_freq(.x, n_cores = n_cores, block_size = block_size, ...)
  # fold the frequencies above 0.5 (which() skips missing values, and only
  # stores the indices of the loci to fold, rather than creating several
  # logical vectors as long as the number of loci)
  to_fold <- which(freq > 0.5)
  freq[to_fold] <- 1 - freq[to_fold]
  freq
}

#' @export
#' @rdname loci_alt_freq
loci_maf.grouped_df <- function(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  type = c("tidy", "list", "matrix"),
  ...
) {
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_maf only works with the genotypes column")
  }

  rlang::check_dots_empty()
  type <- match.arg(type)
  if (is_diploid_only(.x) || is_pseudohaploid(.x)) {
    # matrix of frequencies (loci x groups), without names
    freq_mat <- grouped_alt_freq_dip_pseudo(
      .x,
      as_counts = FALSE,
      n_cores = n_cores,
      block_size = block_size
    )
    # fold the frequencies above 0.5 (which() skips missing values, and only
    # stores the indices of the values to fold)
    to_fold <- which(freq_mat > 0.5)
    freq_mat[to_fold] <- 1 - freq_mat[to_fold]

    groups <- .group_info(.x)
    freq_mat <- format_grouped_output(
      out_mat = freq_mat,
      group_ids = groups$labels,
      loci_names = loci_names(.x),
      type = type,
      group_keys = groups$keys
    )
    return(freq_mat)
  } else {
    # the polyploid case
    stop(
      "loci_maf for polyploid is not directly implemented yet, use",
      "group_map(.x, .f = ~ loci_maf(.x))"
    )
  }
}

# function to estimate frequencies for each group (diploid and pseudohaploid
# data): it returns a plain matrix with one row per locus and one column per
# group of frequencies or, if as_counts = TRUE, with the alternate allele
# counts in the first n_groups columns and the valid allele counts in the
# next n_groups columns
grouped_alt_freq_dip_pseudo <- function(.x,
                                        as_counts = FALSE,
                                        n_cores,
                                        block_size) {
  # the cpp code reads the FBM directly, in blocks of block_size loci, and
  # fills a single matrix (rather than binding the matrices from each
  # block); openMP uses multiple cores in the cpp code. If we want
  # frequencies, it only returns them (freq_only = TRUE)
  grouped_alt_freq_dip_pseudo_cpp(
    BM = .gt_get_fbm(.x),
    rowInd = .gt_fbm_rows(.x),
    colInd = .gt_fbm_cols(.x),
    groupIds = dplyr::group_indices(.x) - 1L,
    ngroups = dplyr::n_groups(.x),
    ploidy = indiv_ploidy(.x),
    ncores = n_cores,
    as_counts = as_counts,
    block_size = block_size,
    freq_only = TRUE
  )
}

# function to estimate frequencies for diploid
loci_alt_freq_dip_pseudo <- function(.x,
                                     as_counts = FALSE,
                                     ploidy,
                                     is_pseudohaploid,
                                     n_cores,
                                     block_size) {
  # get the FBM
  geno_fbm <- .gt_get_fbm(.x)
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # no individuals (e.g. an empty group in group_map()): no valid alleles
  if (length(rows_to_keep) == 0) {
    return(.alt_freq_no_indiv(.x, as_counts))
  }
  # as long as we have more than one individual
  if (length(rows_to_keep) > 1) {
    # the cpp code reads the FBM directly (in blocks of block_size loci, using
    # openMP for parallelisation) and returns either the vector of
    # frequencies, or a matrix with the counts of alternate and valid alleles
    freq <- alt_freq_dip_pseudo_cpp(
      BM = geno_fbm,
      rowInd = rows_to_keep,
      colInd = attr(.x, "loci")$big_index,
      ploidy = ploidy,
      ncores = n_cores,
      as_counts = as_counts,
      block_size = block_size
    )
    return(freq)
  } else {
    # if we have a single individual, create a matrix with two columns
    # counts and ploidy (i.e. number of copies), with no valid alleles (and
    # thus no alternate alleles) for missing genotypes, as in the cpp code
    geno <- geno_fbm[rows_to_keep, attr(.x, "loci")$big_index]
    freq <- cbind(
      n_alt = ifelse(is.na(geno), 0, geno),
      n_valid = ifelse(is.na(geno), 0, ploidy)
    )
    # if this individual is pseudohaploid, get the right counts by dividing by 2
    if (ploidy == 1) {
      freq[, 1] <- freq[, 1] / 2
    }
    # if we are returning frequency, compute it for this individual
    if (!as_counts) {
      freq[, 1] <- freq[, 1] / freq[, 2]
    }
  }
  if (!as_counts) {
    # just return the first column with the frequency
    return(freq[, 1])
  } else {
    return(freq)
  }
}

# The output of loci_alt_freq() when there are no individuals (e.g. an empty
# group in group_map()): no valid alleles at any locus, so the frequencies are
# NA (as for loci without valid alleles), and the counts are zero
.alt_freq_no_indiv <- function(.x, as_counts) {
  n_loci <- length(attr(.x, "loci")$big_index)
  if (as_counts) {
    cbind(n_alt = numeric(n_loci), n_valid = numeric(n_loci))
  } else {
    rep(NA_real_, n_loci)
  }
}

loci_alt_freq_polyploid <- function(.x, as_counts = FALSE, n_cores,
                                    block_size, ...) {
  warning(paste(
    "this function still needs a proper unit test!!!",
    "It assumes alleles are the unit of observation"
  ))
  # get the FBM
  geno_fbm <- .gt_get_fbm(.x)
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # no individuals (e.g. an empty group in group_map()): no valid alleles
  if (length(rows_to_keep) == 0) {
    return(.alt_freq_no_indiv(.x, as_counts))
  }
  # as long as we have more than one individual
  ploidy_by_indiv <- indiv_ploidy(.x)
  if (length(rows_to_keep) > 1) {
    # counts of alternate and valid alleles, accounting for missing values
    # (the cpp code reads the FBM directly, in blocks of block_size loci,
    # using openMP for parallelisation)
    counts <- alt_freq_polyploid_cpp(
      BM = geno_fbm,
      rowInd = rows_to_keep,
      colInd = attr(.x, "loci")$big_index,
      ploidy = ploidy_by_indiv,
      ncores = n_cores,
      as_counts = as_counts,
      block_size = block_size
    )
    # if we only want frequencies, the cpp code returned them directly
    if (!as_counts) {
      return(counts)
    }
  } else {
    # if we have a single individual (ploidy_by_indiv is already restricted
    # to the individuals in .x, so it has length 1); missing genotypes have
    # no valid alleles
    geno <- geno_fbm[rows_to_keep, attr(.x, "loci")$big_index]
    counts <- cbind(
      n_alt = ifelse(is.na(geno), 0, geno),
      n_valid = ifelse(is.na(geno), 0, ploidy_by_indiv)
    )
  }
  if (as_counts) {
    return(counts)
  }
  # now get frequency (unname, as a single locus would keep the column name)
  unname(counts[, 1] / counts[, 2])
}

#' Function to format the matrix for grouped operations
#'
#' This function is used to format the matrix for grouped operations
#'
#' Loci names are only added where they are needed: for "tidy", they become
#' the `loci` column; for "matrix", they become the row names of a
#' [matrix_int_names()] object (so integer names stay integers); "list" has
#' no loci names.
#'
#' @param out_mat a matrix of the output (loci in rows, groups in columns), or
#'   a vector (for a single group or a single locus)
#' @param group_ids a vector of group ids (labels)
#' @param loci_names a vector of loci names
#' @param type a string indicating the type of output
#' @param group_keys optionally, the group keys (a tibble with one row per
#'   group and one column per grouping variable, as [dplyr::group_keys()]);
#'   the "tidy" output then has the values of the grouping variables (keeping
#'   their type, e.g. a factor keeps its levels), in one column per grouping
#'   variable if there are several, rather than the labels in a single
#'   `group` column
#' @returns the data formatted as a tibble, list, or matrix
#' @keywords internal
#' @noRd
format_grouped_output <- function(out_mat, group_ids, loci_names, type,
                                  group_keys = NULL) {
  # use the group labels as names (factors and numbers would otherwise fail
  # as matrix column names, or be used as column positions in tidyselect)
  group_ids <- as.character(group_ids)
  n_loci <- length(loci_names)
  n_groups <- length(group_ids)

  if (type == "tidy") {
    # build the long tibble directly, with one row per locus and group (loci
    # in the outer loop, as given by tidyr::pivot_longer()), rather than
    # casting the matrix to a data.frame and pivoting it, which makes several
    # copies of the data. Once the matrix is transposed, each of its rows
    # (a locus) becomes n_groups consecutive values (this also works if
    # out_mat is a vector, i.e. a single group or a single locus)
    values <- t(out_mat)
    # drop all attributes (dim, dimnames, any class) to get a plain vector
    # (values is not shared, so this happens in place)
    attributes(values) <- NULL
    if (!is.null(group_keys) && ncol(group_keys) > 1) {
      # one column per grouping variable (as dplyr::summarise())
      return(dplyr::bind_cols(
        tibble::tibble(loci = rep(loci_names, each = n_groups)),
        group_keys[rep(seq_len(n_groups), times = n_loci), , drop = FALSE],
        tibble::tibble(value = values)
      ))
    }
    # a single group column, with the values of the grouping variable (keeping
    # its type, e.g. a factor keeps its levels), or the labels if the keys are
    # not given
    if (is.null(group_keys)) {
      group_col <- group_ids
    } else {
      group_col <- group_keys[[1]]
    }
    return(tibble::new_tibble(
      list(
        loci = rep(loci_names, each = n_groups),
        group = rep(group_col, times = n_loci),
        value = values
      ),
      nrow = n_loci * n_groups
    ))
  }

  if (is.vector(out_mat)) {
    # names are added below, so that integer loci names are not cast to
    # characters
    out_mat <- matrix(
      out_mat,
      nrow = n_loci,
      ncol = n_groups
    )
  }
  if (type == "list") {
    # return a list to mimic a group_map
    lapply(seq_len(ncol(out_mat)), function(i) out_mat[, i])
  } else if (type == "matrix") {
    # return a matrix_int_names, so that integer loci names are preserved
    matrix_int_names(out_mat, row_names = loci_names, col_names = group_ids)
  }
}
