# internal functions that make life easier
.gt_fbm_cols <- function(.x) {
  show_loci(.x)$big_index
}

.gt_fbm_rows <- function(.x) {
  vctrs::vec_data(.x$genotypes)
}

.gt_get_fbm <- function(.x) {
  # if this is a gen_tibble
  if (inherits(.x, "gen_tbl")) {
    fbm <- attr(.x$genotypes, "fbm")
  } else {
    fbm <- attr(.x, "fbm")
  }
  if (is.null(fbm)) {
    stop("Missing FBM backing for genotypes; ensure `.x` has attr 'fbm'.")
  }
  fbm
}

# a simple tidier for dist matrices
# tidy.dist is deprecated, and often we have a full matrix rather
# than a dist object
tidy_dist_matrix <- function(mat) {
  if (!inherits(mat, "matrix")) {
    stop("mat should be a matrix")
  }
  xy <- t(utils::combn(colnames(mat), 2))
  colnames(xy) <- c("item1", "item2")
  xy <- xy %>%
    tibble::as_tibble() %>%
    dplyr::mutate(value = mat[xy])
  # add a class to indicate this is a pairwise tibble
  class(xy) <- c("pairwise_tbl", class(xy))
  return(xy)
}

# stop if not diploid
stopifnot_diploid <- function(x) {
  if (inherits(x, "gen_tbl")) {
    x <- x$genotypes
  }
  if (abs(attr(x, "ploidy")) != 2) {
    stop("this function only works on diploid data")
  }
  if (attr(x, "ploidy") == -2) {
    if (min(attr(x, "fbm_ploidy")[vctrs::vec_data(x)]) != 2) {
      stop("this function only works on diploid data")
    }
  }
}

# stop if not diploid or pseudohaploid
stopifnot_dip_pseudo <- function(x) {
  if (inherits(x, "gen_tbl")) {
    x <- x$genotypes
  }
  if (abs(attr(x, "ploidy")) != 2) {
    stop("this function only works on diploid or pseudohaploid data")
  }
}


is_diploid_only <- function(x) {
  if (inherits(x, "gen_tbl")) {
    (attr(x$genotypes, "ploidy") == 2)
  } else {
    (attr(x, "ploidy") == 2)
  }
}

is_pseudohaploid <- function(x) {
  if (inherits(x, "gen_tbl")) {
    (attr(x$genotypes, "ploidy") == -2)
  } else {
    (attr(x, "ploidy") == -2)
  }
}

# stop if the gen_tibble is flagged as containing pseudohaploid data
# (heterozygosity is undefined when only one allele was observed). Trusts
# the ploidy label rather than inspecting which individuals are currently
# selected: dropping pseudohaploid individuals via filter() does not
# refresh the label, so gt_pseudohaploid() must be rerun for this check to
# stop blocking a gen_tibble that no longer actually contains any.
stopifnot_no_pseudohaploid <- function(x) {
  if (is_pseudohaploid(x)) {
    stop(
      "this function does not support pseudohaploid data; ",
      "if you think all pseudohaploid individuals have been ",
      "already removed, rerun gt_pseudohaploid() on your tibble"
    )
  }
}
