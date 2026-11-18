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
  # (there are no pairs with fewer than two individuals)
  n_first <- max(n - 1L, 0L)
  index1 <- rep.int(seq_len(n_first), rev(seq_len(n_first)))
  index2 <- sequence(rev(seq_len(n_first)), from = seq_len(n_first) + 1)
  # (the names of an empty matrix are NULL, which would drop the columns)
  item_names <- if (n == 0) character(0) else colnames(mat)
  xy <- tibble::tibble(
    item1 = item_names[index1],
    item2 = item_names[index2],
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

# evaluate `code` with the temporary FBMs of bigstatsr and bigsnpr in a
# private directory, which is deleted on exit.
#
# bigstatsr's FBM() names its temporary backing files with
# `tempfile(tmpdir = getOption("FBM.dir"))` and appends ".bk". tempfile()
# only checks that the bare name is free, not "<name>.bk", and its random
# suffix has 31 bits, so a new name can match one of the .bk files already in
# the directory (a birthday collision); FBM() then fails with
# "File '<name>.bk' already exists". bigsnpr's snp_clumping() (used by
# snp_autoSVD()) creates a temporary FBM for each chromosome and never
# removes it, so data with many contigs (e.g. the 1341 scaffolds of the
# anolis example) leave thousands of .bk files in tempdir(), and each later
# temporary FBM of the session risks a clash. A fresh, empty directory for
# each call avoids any clash with earlier files, and the leaked files are
# deleted with it. The private directory is made inside the directory set by
# the user for temporary FBMs (`options(FBM.dir)`, `tempdir()` by default),
# e.g. a scratch volume with room for large matrices.
#
# Only use it around code whose value holds no FBM (the FBMs created in the
# private directory are deleted on exit).
# @param code the code to evaluate (lazily, in the caller's frame)
# @returns the value of `code`
# @keywords internal
# @noRd
with_private_fbm_dir <- function(code) {
  local_private_fbm_dir()
  # force the promise now that FBM.dir points to the private directory
  code
}

# the same as with_private_fbm_dir(), for the rest of the calling function:
# the private directory is deleted, and FBM.dir restored, when `envir` exits
# (as withr::local_*() functions). Call it once the function has made any
# FBM that it returns (or only on the paths that return no FBM).
# @param envir the frame whose exit deletes the directory
# @returns the path of the private directory, invisibly
# @keywords internal
# @noRd
local_private_fbm_dir <- function(envir = parent.frame()) {
  # tempfile() checks the bare name, which here is the directory itself, so
  # the directory is new and empty; it is made inside FBM.dir, so that
  # temporary FBMs stay where the user wants them (bigstatsr sets FBM.dir to
  # tempdir() when it is loaded)
  parent_dir <- getOption("FBM.dir", default = tempdir())
  fbm_dir <- tempfile("fbm_tmp_", tmpdir = parent_dir)
  # (an FBM.dir that does not exist, or can not be written to, would
  # otherwise only give an unclear error when the first FBM is made)
  if (!dir.create(fbm_dir, showWarnings = FALSE)) {
    stop(
      "could not create a directory for temporary backing files in '",
      parent_dir, "' (set with options(FBM.dir = ...))",
      call. = FALSE
    )
  }
  old_opts <- options(FBM.dir = fbm_dir)
  # register the clean up in the caller's frame (with add = TRUE, so as not to
  # replace its own on.exit() code); on Windows unlink() can fail on files
  # still mapped, but names drawn in the directory still cannot clash with
  # files elsewhere
  cleanup <- bquote({
    options(.(old_opts))
    unlink(.(fbm_dir), recursive = TRUE)
  })
  do.call(base::on.exit, list(cleanup, add = TRUE), envir = envir)
  invisible(fbm_dir)
}

# Missing genotypes ------------------------------------------------------------
#
# Some analyses (those of bigsnpr and bigstatsr: clumping, randomised and
# partial SVD, pcadapt, matrix products) need complete genotypes. A full pass
# over the genotypes just to check them would cost as much as some of these
# analyses, so the genotypes are checked on a sample, as bigsnpr does
# (stopifnot_no_missing()), and, where an analysis makes a full pass anyway,
# exactly on its results (e.g. scaling_no_missing()). bigsnpr's own sampled
# check is turned off (local_no_bigstatsr_checks()), and all these errors
# come from stop_missing_genotypes(), so that they read the same.

# Stop because the genotypes used by a function have missing values.
# @param what the name of the function that needs complete genotypes (with
#   its brackets), for the message
# @param hint an alternative to imputation, appended to the message
# @param imputed_auto whether the function uses imputed genotypes
#   automatically (otherwise the gen_tibble must be set to use them)
# @keywords internal
# @noRd
stop_missing_genotypes <- function(what, hint = NULL, imputed_auto = TRUE) {
  stop(
    what, " needs complete or imputed genotypes, but some genotypes are ",
    "missing: impute them (e.g. with gt_impute_simple())",
    if (imputed_auto) {
      ", and the imputed genotypes are used automatically"
    } else {
      ", and set the gen_tibble to use them (gt_set_imputed(x, set = TRUE))"
    },
    if (!is.null(hint)) paste0("; ", hint),
    call. = FALSE
  )
}

# Stop if a sample of the genotypes used by x has missing values (after any
# switch to imputed genotypes). Like bigstatsr's assert_noNA(), it reads
# only part of the genotypes, but only the individuals and loci that x uses
# (the backing matrix can also hold individuals and loci, with missing
# genotypes, that x no longer uses): all the individuals, at evenly spaced
# loci (no random numbers, so that the user's seed is left alone), for about
# 1e6 genotypes and at least 100 loci. Missing genotypes at a few loci can
# escape it, so the help pages ask for all of them to be imputed (or the
# loci removed) first.
# @param x a gen_tibble, or its genotypes column
# @param what,hint as in stop_missing_genotypes()
# @param loci the indices of the loci of x that the function uses (all of
#   them by default)
# @param n_cores the number of cores
# @keywords internal
# @noRd
stopifnot_no_missing <- function(x, what, hint = NULL, loci = NULL,
                                 n_cores = 1) {
  if (!inherits(x, "vctrs_bigSNP")) {
    stopifnot_gen_tibble(x)
    x <- x$genotypes
  }
  # the rows (individuals) and columns (loci) of the backing matrix used by x
  rows <- vctrs::vec_data(x)
  cols <- attr(x, "loci")$big_index
  if (!is.null(loci)) {
    cols <- cols[loci]
  }
  if (length(rows) == 0 || length(cols) == 0) {
    return(invisible(TRUE))
  }
  # a sample of evenly spaced loci, with all the individuals
  n_sample <- min(length(cols), max(100, ceiling(1e6 / length(rows))))
  cols <- cols[unique(round(seq(1, length(cols), length.out = n_sample)))]
  # count the missing genotypes at each sampled locus
  n_na <- missingness_cpp(
    BM = .gt_get_fbm(x), rowInd = rows, colInd = cols, ncores = n_cores
  )
  if (any(n_na > 0)) {
    stop_missing_genotypes(what, hint = hint)
  }
  invisible(TRUE)
}

# Turn off the argument checks of bigstatsr and bigsnpr (their option
# "bigstatsr.check.args") until `envir` exits. Their check for missing
# values (assert_noNA()) samples the whole backing matrix, not the rows and
# columns that the analysis uses, and so gives false errors on a gen_tibble
# that no longer uses some individuals or loci with missing genotypes; the
# genotypes are checked with stopifnot_no_missing() instead. (The other
# checks are of indices and numbers of cores, which tidypopgen builds.)
# @param envir the frame whose exit restores the option
# @keywords internal
# @noRd
local_no_bigstatsr_checks <- function(envir = parent.frame()) {
  old_opts <- options(bigstatsr.check.args = FALSE)
  # restore the option when envir exits (add = TRUE keeps its own on.exit()
  # code)
  do.call(base::on.exit, list(bquote(options(.(old_opts))), add = TRUE),
    envir = envir
  )
  invisible(old_opts)
}

# A scaling function for the SVDs of bigstatsr and bigsnpr (big_SVD(),
# big_randomSVD() and snp_autoSVD()) that first turns off the argument checks
# of bigstatsr (see local_no_bigstatsr_checks()). With n_cores > 1,
# big_randomSVD() computes the scaling and the matrix products in the workers
# of a PSOCK cluster, which do not inherit the options of the session, and
# the scaling is the first thing each worker runs; there, bigstatsr's
# sampled check for missing values in the products could fail on genotypes
# that the gen_tibble does not use, and an error in a worker leaves
# big_randomSVD() waiting forever. In the session itself, the option is
# already off (and restored by local_no_bigstatsr_checks()).
# @param fun_scaling a scaling function, e.g. bigsnpr::snp_scaleBinom()
# @returns a scaling function with the same arguments
# @keywords internal
# @noRd
scaling_without_checks <- function(fun_scaling) {
  scaling <- function(X, ind.row, ind.col) { # nolint: object_name_linter
    # (in a worker, the option stays off until the cluster is stopped)
    options(bigstatsr.check.args = FALSE)
    fun_scaling(X, ind.row = ind.row, ind.col = ind.col)
  }
  # the function is sent to the workers: keep its environment small, without
  # the namespace of tidypopgen (which the workers would have to load)
  environment(scaling) <- list2env(
    list(fun_scaling = fun_scaling),
    parent = baseenv()
  )
  scaling
}
