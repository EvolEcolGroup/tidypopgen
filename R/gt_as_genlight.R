#' Convert a `gen_tibble` to a `genlight` object from `adegenet`
#'
#' This function converts a `gen_tibble` to a `genlight` object from `adegenet`
#'
#' The chromosome and position of each locus are stored in the `genlight`
#' object, and the alleles are stored as "ref/alt" in `loc.all` (`genlight`
#' counts the second allele, i.e. the alternate allele). As `genlight` only
#' supports single letter alleles, `loc.all` is not set (with a warning) if any
#' locus has a missing or multi-character allele.
#'
#' The genotypes are read in chunks of loci, so that the full matrix of
#' genotypes is never loaded in memory; the result is the same as building the
#' `genlight` from that matrix.
#'
#' @param x a [`gen_tibble`], with population coded as 'population'
#' @param chunk_size the number of loci processed at a time (rounded up to a
#'   multiple of 8). If left to NULL, chunks of about 10 million genotypes
#'   are used.
#' @returns a `genlight` object
#' @export
#' @examplesIf rlang::is_installed("adegenet")
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Convert to genlight
#' gt_genlight <- example_gt %>% gt_as_genlight()
#'
#' # Check object class
#' class(gt_genlight)
#'
gt_as_genlight <- function(x, chunk_size = NULL) {
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop(
      "to use this function, first install package 'adegenet' with\n",
      "install.packages('adegenet')"
    )
  }
  stopifnot_diploid(x) # currently we only support diploid data
  loci <- show_loci(x)
  # genlight counts the second allele in loc.all, i.e. the alternate allele
  loc_all <- paste0(loci$allele_ref, "/", loci$allele_alt)
  # adegenet only accepts single letter alleles (e.g. "a/t") for all loci
  if (!all(grepl("^[[:alpha:]]{1}/[[:alpha:]]{1}$", loc_all))) {
    warning(
      "some loci have missing or multi-character alleles, which genlight ",
      "does not support; the alleles (loc.all) were not stored"
    )
    loc_all <- NULL
  }
  test_genlight <- methods::new( # nolint
    "genlight", # nolint
    gen = gt_as_snpbin_list(x, chunk_size = chunk_size),
    ploidy = 2, # TODO update this when ploidy is implemented
    ind.names = x$id,
    pop = x$population,
    loc.names = loci_names(x),
    loc.all = loc_all,
    chromosome = loci$chromosome,
    position = loci$position
  )
}

# the genotypes of each individual of a diploid gen_tibble as a list of
# `SNPbin` objects (from `adegenet`), reading the genotypes in chunks of loci.
# The objects are the same as those that new("SNPbin") makes from the vector
# of dosages of each individual (as new("genlight") does from a matrix): the
# bits of loci with a dosage of 2, then those with a dosage of 1 or 2 (only
# the latter if the individual has no dosage of 2), with missing genotypes as
# 0 and their positions in NA.posi, and the largest dosage as ploidy (at least
# 1, or NA if all the genotypes of the individual are missing).
gt_as_snpbin_list <- function(x, chunk_size = NULL) {
  n_indiv <- nrow(x)
  n_loci <- count_loci(x)
  # by default, chunks of about 10 million genotypes
  if (is.null(chunk_size)) {
    chunk_size <- 1e7 %/% n_indiv
  } else if (!is.numeric(chunk_size) || length(chunk_size) != 1 ||
               is.na(chunk_size) || chunk_size < 1) {
    stop("'chunk_size' should be a positive number")
  }
  # chunks hold a multiple of 8 loci, so that their bits fill whole bytes (a
  # chunk larger than the number of loci, e.g. Inf, holds all the loci)
  chunk_size <- min(chunk_size, n_loci)
  chunk_size <- max(8, ceiling(chunk_size / 8) * 8)
  n_chunks <- ceiling(n_loci / chunk_size)

  # bits of each individual (one column each), packed 8 loci per byte: loci
  # with a dosage of 2, and loci with a dosage of 1 or 2
  n_bytes <- ceiling(n_loci / 8)
  bits_two <- matrix(as.raw(0), nrow = n_bytes, ncol = n_indiv)
  bits_any <- matrix(as.raw(0), nrow = n_bytes, ncol = n_indiv)
  # whether each individual has a dosage of 2, and its number of genotypes
  # that are not missing
  has_two <- rep(FALSE, n_indiv)
  n_called <- rep(0, n_indiv)
  # positions of the missing genotypes of each individual, by chunk
  na_posi_chunks <- vector("list", n_chunks)

  for (chunk_i in seq_len(n_chunks)) {
    chunk_start <- (chunk_i - 1) * chunk_size + 1
    this_chunk <- chunk_start:min(chunk_start + chunk_size - 1, n_loci)
    # dosages of this chunk, with loci in rows, as integers (as in
    # new("SNPbin"), which truncates them)
    geno <- t(show_genotypes(x, loci_indices = this_chunk))
    storage.mode(geno) <- "integer"
    is_missing <- is.na(geno)
    n_called <- n_called + length(this_chunk) - colSums(is_missing)
    # positions of the missing genotypes, for each individual (the matrix
    # has one column per individual); the factor of individuals is built
    # from its codes, as factor() would turn each of them into a string
    missing_index <- which(is_missing) - 1L
    missing_indiv <- structure(
      missing_index %/% length(this_chunk) + 1L,
      levels = as.character(seq_len(n_indiv)),
      class = "factor"
    )
    na_posi_chunks[[chunk_i]] <- split(
      missing_index %% length(this_chunk) + as.integer(chunk_start),
      missing_indiv
    )
    # missing genotypes are not set in either set of bits
    geno[is_missing] <- 0L
    is_two <- geno == 2L
    is_any <- geno >= 1L
    has_two <- has_two | colSums(is_two) > 0
    # pad the last chunk with loci that are not set, to fill its last byte
    n_pad <- (8 - length(this_chunk) %% 8) %% 8
    if (n_pad > 0) {
      is_two <- rbind(is_two, matrix(FALSE, nrow = n_pad, ncol = n_indiv))
      is_any <- rbind(is_any, matrix(FALSE, nrow = n_pad, ncol = n_indiv))
    }
    # pack the bits of each individual (a column) into bytes
    chunk_bytes <- (chunk_start - 1) / 8 + seq_len(nrow(is_two) / 8)
    bits_two[chunk_bytes, ] <- packBits(is_two)
    bits_any[chunk_bytes, ] <- packBits(is_any)
  }

  # create a SNPbin for each individual
  lapply(seq_len(n_indiv), function(i) {
    na_posi <- unlist(lapply(na_posi_chunks, `[[`, i), use.names = FALSE)
    if (n_called[i] == 0) {
      # all genotypes are missing
      snp <- list(bits_any[, i])
      ploidy <- NA_integer_
    } else if (has_two[i]) {
      snp <- list(bits_two[, i], bits_any[, i])
      ploidy <- 2L
    } else {
      snp <- list(bits_any[, i])
      ploidy <- 1L
    }
    methods::new(
      "SNPbin",
      snp = snp,
      n.loc = n_loci,
      NA.posi = na_posi,
      ploidy = ploidy
    )
  })
}
