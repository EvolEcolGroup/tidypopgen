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
  n <- ncol(mat)
  # indices of all pairs (i < j), in the same order as utils::combn(n, 2),
  # i.e. (1,2), (1,3), ..., (1,n), (2,3), ...; using integer indices avoids
  # building a character matrix of names and matching names to index mat
  index1 <- rep.int(seq_len(n - 1), rev(seq_len(n - 1)))
  index2 <- sequence(rev(seq_len(n - 1)), from = seq_len(n - 1) + 1)
  xy <- tibble::tibble(
    item1 = colnames(mat)[index1],
    item2 = colnames(mat)[index2],
    value = mat[cbind(index1, index2)]
  )
  # add a class to indicate this is a pairwise tibble
  class(xy) <- c("pairwise_tbl", class(xy))
  return(xy)
}

# stop if not diploid
stopifnot_diploid <- function(x) {
  if (inherits(x, "tbl_df")) {
    stopifnot_gen_tibble(x)
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
  if (inherits(x, "tbl_df")) {
    stopifnot_gen_tibble(x)
    x <- x$genotypes
  }
  if (abs(attr(x, "ploidy")) != 2) {
    stop("this function only works on diploid or pseudohaploid data")
  }
}


is_diploid_only <- function(x) {
  if (inherits(x, "tbl_df")) {
    stopifnot_gen_tibble(x)
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

# refresh a gen_tibble's top-level ploidy label after some individuals have
# been dropped (e.g. by filter()), which does not update it automatically.
# Only two cases can safely be inferred purely from which individuals
# remain:
# - pseudohaploid (-2): reset to diploid (2) if no pseudohaploid
#   individuals remain (pseudohaploid is a data-quality flag, not a
#   ploidy value, so it is only ever reset to 2, never to whatever
#   fbm_ploidy the remaining individuals happen to have)
# - mixed ploidy (0): reset to that single ploidy if the remaining
#   individuals now all happen to share the same one
# any other ploidy value is a definite single ploidy (diploid or
# polyploid) that dropping individuals cannot make stale, so it is left
# untouched.
.gt_refresh_ploidy <- function(out) {
  remaining_ploidy <-
    attr(out$genotypes, "fbm_ploidy")[vctrs::vec_data(out$genotypes)]
  if (length(remaining_ploidy) == 0) {
    return(out)
  }
  if (is_pseudohaploid(out)) {
    if (min(remaining_ploidy) == 2) {
      attr(out$genotypes, "ploidy") <- 2
    }
  } else if (attr(out$genotypes, "ploidy") == 0) {
    if (length(unique(remaining_ploidy)) == 1) {
      attr(out$genotypes, "ploidy") <- remaining_ploidy[1]
    }
  }
  out
}

# stop if the gen_tibble is flagged as containing pseudohaploid data
# (heterozygosity is undefined when only one allele was observed). Trusts
# the ploidy label rather than inspecting which individuals are currently
# selected: filter() refreshes the label where it safely can (see
# .gt_refresh_ploidy()), but e.g. dropping only some pseudohaploid
# individuals, or subsetting some other way, will not, so rerun
# gt_pseudohaploid() if this keeps erroring unexpectedly.
stopifnot_no_pseudohaploid <- function(x) {
  if (is_pseudohaploid(x)) {
    stop(
      "this function does not support pseudohaploid data; ",
      "if you think all pseudohaploid individuals have been ",
      "already removed, rerun gt_pseudohaploid() on your tibble"
    )
  }
}
