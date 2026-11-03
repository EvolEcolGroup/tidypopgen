#' Compute population specific Fst
#'
#' This function computes population specific Fst, using the approach in Weir
#' and Goudet 2017 (as computed by [hierfstat::fst.dosage()]).
#' @references Weir, BS and Goudet J (2017) A Unified Characterization of
#'   Population Structure and Relatedness. Genetics (2017) 206:2085
#' @param .x a grouped [`gen_tibble`] (as obtained by using
#'   [dplyr::group_by()]). If grouped by more than one variable, groups are
#'   labelled by the values of the grouping variables pasted with "." (e.g.
#'   `pop1.F`).
#' @param include_global boolean determining whether, besides the population
#'   specific Fst, a global Fst should be appended. Note that this will return a
#'   vector of n populations plus 1 (the global value)
#' @param allele_sharing_mat optional, the matrix of Allele Sharing returned by
#'   [pairwise_allele_sharing()] with `as_matrix=TRUE`. As a number of
#'   statistics can be derived from the Allele Sharing matrix,
#' @returns a vector of population specific Fst (plus the global value if
#'   `include_global=TRUE`)
#' @export
#' @seealso [hierfstat::fst.dosage()]
#' @examples
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # Compute FIS using Nei87
#' example_gt %>% pop_fst()
#'
#' # To include the global Fst, set include_global = TRUE
#' example_gt %>% pop_fst(include_global = TRUE)
#'
#' # To calculate from a pre-computed allele sharing matrix:
#' allele_sharing_mat <- pairwise_allele_sharing(example_gt, as_matrix = TRUE)
#' example_gt %>% pop_fst(allele_sharing_mat = allele_sharing_mat)
#'
pop_fst <- function(.x, include_global = FALSE, allele_sharing_mat = NULL) {
  stopifnot_diploid(.x)
  if (!inherits(.x, "grouped_df")) {
    stop(".x should be a grouped gen_tibble")
  }
  if (is.null(allele_sharing_mat)) {
    allele_sharing_mat <- pairwise_allele_sharing(.x, as_matrix = TRUE)
  }

  Mij <- allele_sharing_mat
  diag(Mij) <- NA
  groups <- .group_info(.x)
  # the group of each individual
  pop <- groups$ids + 1L
  # empty groups (which only exist if the user asked for them) do not
  # contribute to the mean allele sharing between groups, and get NaN
  pop_levels <- which(groups$sizes > 0)
  n_pop <- length(pop_levels)
  if (n_pop < 2) {
    stop("pop_fst requires at least two groups")
  }
  wil <- lapply(pop_levels, function(z) which(pop == z))
  Fsts <- unlist(lapply(
    wil,
    function(pop_levels) {
      mean(Mij[pop_levels, pop_levels], na.rm = TRUE)
    }
  ))
  Mb <- 0
  mMij <- matrix(numeric(n_pop^2), ncol = n_pop)
  for (i in 2:n_pop) {
    p1 <- wil[[i]]
    for (j in 1:(i - 1)) {
      p2 <- wil[[j]]
      mMij[i, j] <- mMij[j, i] <- mean(Mij[p1, p2], na.rm = TRUE)
      Mb <- Mb + mMij[i, j]
    }
  }
  diag(mMij) <- Fsts

  Mb <- Mb * 2 / (n_pop * (n_pop - 1))

  # pop specific fsts
  fst_non_empty <- (Fsts - Mb) / (1 - Mb)
  fst_by_pop <- rep(NaN, groups$n)
  fst_by_pop[pop_levels] <- fst_non_empty

  if (include_global) {
    fst_by_pop <- c(fst_by_pop, mean(fst_non_empty, na.rm = TRUE))
    names(fst_by_pop) <- c(groups$labels, "global")
  } else {
    names(fst_by_pop) <- groups$labels
  }
  return(fst_by_pop)
}
