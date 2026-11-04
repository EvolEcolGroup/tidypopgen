# A reference implementation of the site frequency spectra of easySFS, used
# to test pop_sfs(), pairwise_pop_sfs(), nwise_pop_sfs() and
# pop_sfs_preview(). It is a literal port of easySFS's make_datadict() and of
# the functions of dadi's Spectrum class bundled with easySFS
# (count_data_dict(), _from_count_dict(), _cached_projection(), fold() and
# S()), from https://github.com/isaacovercast/easySFS (commit 87d2996).
# Spectra are returned as R arrays with one dimension per population (so the
# first population changes fastest in memory), with a "mask" attribute.

# make_datadict(): the number of called and alternate allele copies of each
# locus in each population. alt_copies is a matrix (individuals x loci) of the
# number of alternate allele copies of each individual (NA if missing), and
# ploidy the number of allele copies of each individual (one for
# pseudohaploids). Returns a list (one element per population in pops) of
# matrices with columns called and derived (one row per locus).
ref_sfs_calls <- function(alt_copies, ploidy, population, pops) {
  lapply(pops, function(pop) {
    rows <- population == pop
    alt <- alt_copies[rows, , drop = FALSE]
    cbind(
      called = colSums((!is.na(alt)) * ploidy[rows]),
      derived = colSums(alt, na.rm = TRUE)
    )
  })
}

# count_data_dict(): the number of loci with each configuration of called and
# derived copies in each population (a data.frame with columns called_1, ...,
# called_K, derived_1, ..., derived_K and count)
ref_sfs_count_dict <- function(calls) {
  configs <- do.call(cbind, c(
    lapply(calls, function(x) x[, "called"]),
    lapply(calls, function(x) x[, "derived"])
  ))
  k <- length(calls)
  colnames(configs) <- c(paste0("called_", seq_len(k)), paste0(
    "derived_",
    seq_len(k)
  ))
  configs <- as.data.frame(configs)
  configs$count <- 1
  stats::aggregate(count ~ ., data = configs, FUN = sum)
}

# _lncomb() and _cached_projection(): the probability of each number of
# derived copies (0 to proj_to) when projecting proj_from copies, hits of
# which are derived, down to proj_to copies
ref_sfs_lncomb <- function(n, k) {
  lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
}
ref_sfs_projection <- function(proj_to, proj_from, hits) {
  contrib <- rep(0, proj_to + 1)
  if (proj_from < proj_to) {
    return(contrib)
  }
  proj_hits <- 0:proj_to
  # scipy's gammaln() is Inf for negative integers, so that the terms with
  # more derived copies than available (hits - proj_hits < 0), or more
  # ancestral copies than available, are 0
  valid <- (hits - proj_hits) >= 0 &
    (proj_from - proj_to) - (hits - proj_hits) >= 0
  proj_hits <- proj_hits[valid]
  lncontrib <- ref_sfs_lncomb(proj_to, proj_hits) +
    ref_sfs_lncomb(proj_from - proj_to, hits - proj_hits) -
    ref_sfs_lncomb(proj_from, hits)
  contrib[valid] <- exp(lncontrib)
  contrib
}

# reverse_array(): reverse an array along all its axes
ref_sfs_reverse <- function(x) {
  idx <- lapply(dim(x), function(d) rev(seq_len(d)))
  do.call(`[`, c(list(x), idx, list(drop = FALSE)))
}

# fold(): fold a spectrum on the total number of derived copies
ref_sfs_fold <- function(fs) {
  mask <- attr(fs, "mask")
  dims <- dim(fs)
  total_samples <- sum(dims - 1)
  # the total number of derived copies of each cell
  idx <- as.matrix(expand.grid(lapply(dims, function(d) seq_len(d) - 1)))
  total_per_entry <- array(rowSums(idx), dim = dims)
  where_folded_out <- total_per_entry > floor(total_samples / 2)
  original_mask <- mask
  final_mask <- original_mask | ref_sfs_reverse(original_mask)
  data <- array(as.vector(fs), dim = dims)
  reversed <- ref_sfs_reverse(array(
    ifelse(where_folded_out, data, 0),
    dim = dims
  ))
  folded <- data + reversed
  folded[where_folded_out] <- 0
  where_ambiguous <- total_per_entry == total_samples / 2
  ambiguous <- array(ifelse(where_ambiguous, data, 0), dim = dims)
  folded <- folded - 0.5 * ambiguous + 0.5 * ref_sfs_reverse(ambiguous)
  final_mask <- final_mask | where_folded_out
  structure(array(folded, dim = dims), mask = final_mask)
}

# Spectrum.from_data_dict(): the spectrum of the populations in pops,
# projected to proj copies, folded unless polarized = TRUE (the REF allele is
# the outgroup allele in easySFS, so all loci are polarised)
ref_sfs <- function(alt_copies, ploidy, population, pops, proj,
                    polarized = FALSE) {
  calls <- ref_sfs_calls(alt_copies, ploidy, population, pops)
  configs <- ref_sfs_count_dict(calls)
  k <- length(pops)
  fs <- array(0, dim = proj + 1)
  for (i in seq_len(nrow(configs))) {
    contribs <- lapply(seq_len(k), function(i_pop) {
      ref_sfs_projection(
        proj[i_pop],
        configs[i, paste0("called_", i_pop)],
        configs[i, paste0("derived_", i_pop)]
      )
    })
    fs_proj <- array(Reduce(outer, contribs), dim = proj + 1)
    fs <- fs + configs$count[i] * fs_proj
  }
  # mask the corners
  mask <- array(FALSE, dim = proj + 1)
  mask[1] <- TRUE
  mask[length(mask)] <- TRUE
  fs <- structure(fs, mask = mask)
  if (polarized) {
    fs
  } else {
    ref_sfs_fold(fs)
  }
}

# Spectrum.S(): the number of segregating sites (the sum of the cells that
# are neither masked nor corners)
ref_sfs_seg_sites <- function(fs) {
  mask <- attr(fs, "mask")
  mask[1] <- TRUE
  mask[length(mask)] <- TRUE
  sum(fs[!mask])
}

# simulate a grouped gen_tibble with ragged missingness; returns a list with
# the gen_tibble (gt), the matrix of alternate allele copies (alt_copies),
# the ploidy of each individual and its population. modify is an optional
# function of the genotypes and populations, returning modified genotypes
# (e.g. to add missing data)
ref_sfs_simulate <- function(pop_sizes, n_loci, ploidy = 2,
                             max_missing = 0.5, seed = 1, modify = NULL) {
  withr::local_seed(seed)
  population <- rep(names(pop_sizes), pop_sizes)
  n_indiv <- length(population)
  ploidy <- rep_len(ploidy, n_indiv)
  anc_freq <- stats::rbeta(n_loci, 0.3, 0.3)
  pop_freq <- sapply(names(pop_sizes), function(pop) {
    stats::rbeta(n_loci, anc_freq * 20 + 0.01, (1 - anc_freq) * 20 + 0.01)
  })
  genotypes <- matrix(NA_real_, nrow = n_indiv, ncol = n_loci)
  for (i_ind in seq_len(n_indiv)) {
    genotypes[i_ind, ] <- stats::rbinom(
      n_loci,
      size = ploidy[i_ind],
      prob = pop_freq[, population[i_ind]]
    )
  }
  missing_rate <- stats::runif(n_loci, 0, max_missing)
  for (i_locus in seq_len(n_loci)) {
    genotypes[stats::runif(n_indiv) < missing_rate[i_locus], i_locus] <- NA
  }
  if (!is.null(modify)) {
    genotypes <- modify(genotypes, population)
  }
  loci <- data.frame(
    name = paste0("rs", seq_len(n_loci)),
    chromosome = 1,
    position = seq_len(n_loci),
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
  indiv_meta <- data.frame(
    id = paste0("ind", seq_len(n_indiv)),
    population = population
  )
  args <- list(
    x = genotypes,
    loci = loci,
    indiv_meta = indiv_meta,
    backingfile = tempfile(),
    quiet = TRUE
  )
  if (any(ploidy != 2)) {
    args$ploidy <- ploidy
  }
  gt <- do.call(gen_tibble, args) %>% dplyr::group_by(population)
  list(
    gt = gt,
    alt_copies = genotypes,
    ploidy = ploidy,
    population = population
  )
}

# the values (without attributes) of an sfs object or a reference spectrum
sfs_values <- function(fs) {
  array(as.vector(fs), dim = dim(fs))
}

# read a spectrum written by dadi's Spectrum.to_file() (as used by easySFS):
# a line with the dimensions, "folded" or "unfolded" and the population
# labels, a line with the cells (the last population changing fastest), and
# a line with the mask
read_dadi_sfs <- function(file) {
  lines <- readLines(file)
  lines <- lines[!startsWith(lines, "#")]
  header <- strsplit(trimws(lines[1]), " +")[[1]]
  is_dim <- grepl("^[0-9]+$", header)
  list(
    dims = as.integer(header[is_dim]),
    folded = header[!is_dim][1] == "folded",
    pop_ids = gsub("\"", "", header[!is_dim][-1]),
    data = as.numeric(strsplit(trimws(lines[2]), " +")[[1]]),
    mask = as.logical(as.integer(strsplit(trimws(lines[3]), " +")[[1]]))
  )
}

# the text of a file, with its numbers replaced by "#", and its numbers (in
# the order they appear), to compare files whose numbers may differ in their
# last digits
sfs_file_parts <- function(file) {
  lines <- readLines(file)
  num_re <- "-?[0-9]+(\\.[0-9]+)?(e[-+]?[0-9]+)?"
  list(
    text = gsub(num_re, "#", lines),
    numbers = as.numeric(unlist(regmatches(lines, gregexpr(num_re, lines))))
  )
}
