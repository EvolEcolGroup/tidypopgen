#' Simulate a toy dataset of unlinked loci
#'
#' Simulates the genotypes of a few populations at unlinked biallelic loci
#' under the Balding-Nichols model, and returns them as a `gen_tibble`. This is
#' a **toy simulator** for quick checks of methods on unlinked loci (e.g. that
#' a PCA separates populations, or how a statistic behaves with missing
#' genotypes); it is not meant for inference.
#'
#' Each locus has an ancestral allele frequency \eqn{p} drawn uniformly in
#' `freq_range`. The frequency in population \eqn{k} is drawn from a beta
#' distribution with mean \eqn{p} and variance \eqn{F_k p (1-p)} (the
#' Balding-Nichols model, i.e. independent drift of each population from a
#' common ancestor, a star-shaped tree), and the genotypes of its individuals
#' are binomial draws of two alleles at that frequency (Hardy-Weinberg
#' proportions). Genotypes are then set to missing at random with
#' probability `missing`.
#'
#' The model has no linkage, mutation, recombination, demographic history
#' or tree of populations: for realistic data, simulate with the coalescent
#' (e.g. msprime, which can be run from R with the slendr package).
#'
#' All loci are placed on chromosome 1, 1,000 bp apart, but they are
#' independent. Some loci can be monomorphic in the sample. The genotypes are
#' generated in memory before being written to the backing file, so keep the
#' datasets small. Use [set.seed()] beforehand for reproducible datasets.
#'
#' @param n_indiv the number of individuals per population: one value, or one
#'   per population.
#' @param n_loci the number of loci.
#' @param fst the divergence (\eqn{F_{ST}}) of each population from the
#'   ancestral population: one value, or one per population. The number of
#'   populations is the longer of `n_indiv` and `fst`.
#' @param missing the probability that a genotype is missing: one value, or
#'   one per individual (e.g. to give individuals different missingness).
#' @param freq_range the range of the (uniform) distribution of the ancestral
#'   allele frequencies.
#' @param backingfile the path, without extension, of the backing files; by
#'   default a temporary file.
#' @returns a `gen_tibble` with columns `id`, `population` and `genotypes`.
#' @references Balding D.J. and Nichols R.A. (1995) A method for quantifying
#'   differentiation between populations at multi-allelic loci and its
#'   implications for investigating identity and paternity. Genetica 96: 3-12.
#'
#'   Baumdicker F. et al. (2022) Efficient ancestry and mutation simulation
#'   with msprime 1.0. Genetics 220: iyab229.
#'
#'   Petr M., Haller B.C., Ralph P.L. and Racimo F. (2023) slendr: a framework
#'   for spatio-temporal population genomic simulations on geographic
#'   landscapes. Peer Community Journal 3: e121.
#' @export
#' @examples
#' set.seed(1)
#' # three populations of 20 individuals, 500 loci
#' sim_gt <- gt_simulate_bn(
#'   n_indiv = 20, n_loci = 500, fst = c(0.02, 0.05, 0.1)
#' )
#' sim_gt
#' # with missing genotypes, increasing from 0 to 90% across individuals
#' sim_missing <- gt_simulate_bn(
#'   n_indiv = 20, n_loci = 500, fst = 0.05,
#'   missing = rep(seq(0, 0.9, length.out = 20), 3)
#' )
gt_simulate_bn <- function(n_indiv = 20, n_loci = 1000, fst = 0.05,
                           missing = 0, freq_range = c(0.05, 0.95),
                           backingfile = NULL) {
  # the number of populations, and the values for each of them
  n_pops <- max(length(n_indiv), length(fst))
  if (!length(n_indiv) %in% c(1, n_pops) || !length(fst) %in% c(1, n_pops)) {
    stop("'n_indiv' and 'fst' should have one value, or one per population")
  }
  n_indiv <- rep_len(n_indiv, n_pops)
  fst <- rep_len(fst, n_pops)
  if (any(n_indiv < 1) || any(n_indiv != round(n_indiv))) {
    stop("'n_indiv' should be positive integers")
  }
  if (length(n_loci) != 1 || n_loci < 1 || n_loci != round(n_loci)) {
    stop("'n_loci' should be a positive integer")
  }
  if (any(fst < 0) || any(fst >= 1)) {
    stop("'fst' should be in [0, 1)")
  }
  bad_range <- length(freq_range) != 2 || freq_range[1] < 0 ||
    freq_range[2] > 1 || freq_range[1] > freq_range[2]
  if (bad_range) {
    stop("'freq_range' should be two frequencies, the lower first")
  }
  n_total <- sum(n_indiv)
  bad_missing <- !length(missing) %in% c(1, n_total) || any(missing < 0) ||
    any(missing > 1)
  if (bad_missing) {
    stop("'missing' should be probabilities: one, or one per individual")
  }

  # ancestral allele frequencies
  p_anc <- stats::runif(n_loci, freq_range[1], freq_range[2])
  # population frequencies, drifted independently from the ancestor
  pop_freq <- lapply(fst, function(f) {
    if (f == 0) {
      return(p_anc)
    }
    stats::rbeta(n_loci, p_anc * (1 - f) / f, (1 - p_anc) * (1 - f) / f)
  })
  # genotypes (alternate allele counts), individuals by loci
  genotypes <- do.call(rbind, lapply(seq_len(n_pops), function(k) {
    # each individual of population k, at each locus, at its frequency
    freq_k <- rep(pop_freq[[k]], each = n_indiv[k])
    matrix(stats::rbinom(n_indiv[k] * n_loci, 2, freq_k), nrow = n_indiv[k])
  }))
  # missing genotypes, with a probability per individual
  missing <- rep_len(missing, n_total)
  if (any(missing > 0)) {
    is_missing <- matrix(stats::runif(n_total * n_loci), n_total) < missing
    genotypes[is_missing] <- NA
  }

  indiv_meta <- data.frame(
    id = paste0("ind", seq_len(n_total)),
    population = rep(paste0("pop", seq_len(n_pops)), times = n_indiv)
  )
  loci <- data.frame(
    name = paste0("snp", seq_len(n_loci)),
    chromosome = "1",
    position = seq_len(n_loci) * 1000L,
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "G"
  )
  if (is.null(backingfile)) backingfile <- tempfile("gt_simulate_bn_")
  gen_tibble(genotypes,
    indiv_meta = indiv_meta, loci = loci,
    backingfile = backingfile, quiet = TRUE
  )
}
