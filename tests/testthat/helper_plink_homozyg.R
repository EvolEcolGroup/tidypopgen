# A slow R transcription of PLINK 1.9 --homozyg (calc_homozyg() and
# roh_update() in plink_homozyg.c, without the extend modifier), for one
# individual, used to test windows_indiv_roh().
#
# geno: the dosages of the individual (0, 1, 2 or NA); chrom, pos: the
# chromosome and position of each locus (sorted). The other arguments are
# those of PLINK, in its units (kb). Returns a data frame with the first and
# last locus (indices into geno), the number of loci and the numbers of
# heterozygous and missing genotypes of each run.
plink_homozyg_r <- function(geno, chrom, pos, window_snp = 50,
                            window_het = 1, window_missing = 5,
                            window_threshold = 0.05, min_snp = 100,
                            min_kb = 1000, max_kb_per_snp = 50,
                            max_gap_kb = 1000, max_het = Inf,
                            subtract_1_from_lengths = FALSE) {
  eps <- 2^-30
  min_bases <- floor(min_kb * 1000 + eps)
  max_bases_per_snp <- max_kb_per_snp * 1000 + eps
  max_gap <- floor(max_gap_kb * 1000 + eps)
  new_lengths <- if (subtract_1_from_lengths) 0 else 1
  w <- window_snp
  runs <- list()
  for (ch in unique(chrom)) {
    idx <- which(chrom == ch)
    n <- length(idx)
    if (n < w) next
    g <- geno[idx]
    p <- pos[idx]
    is_het <- !is.na(g) & g == 1
    is_miss <- is.na(g)
    # the hit of the window that starts at each locus (0-based c = 0..n-w)
    win_hit <- vapply(seq_len(n - w + 1), function(s) {
      win <- s:(s + w - 1)
      sum(is_het[win]) <= window_het && sum(is_miss[win]) <= window_missing
    }, logical(1))
    open <- FALSE
    start <- 0
    n_het <- 0
    n_miss <- 0
    # c is 0-based, as in PLINK; c = n closes the last run
    for (c in 0:n) {
      if (c < n) {
        # windows that start in [c - w + 1, min(c, n - w)]
        first_win <- max(0, c - w + 1)
        last_win <- min(c, n - w)
        swhit <- sum(win_hit[(first_win:last_win) + 1])
        n_thr <- if (c <= n - w) min(c + 1, w) else n - c
        swhit_min <- as.integer(n_thr * window_threshold + 1 - eps)
        hit <- swhit >= swhit_min
        call_het <- is_het[c + 1]
        forced_end <- c > 0 && (p[c + 1] - p[c]) > max_gap
      } else {
        hit <- FALSE
        call_het <- FALSE
        forced_end <- FALSE
      }
      if (open && (!hit || (call_het && n_het == max_het) || forced_end)) {
        len <- c - start
        base_len <- p[c] + new_lengths - p[start + 1]
        keep <- len >= min_snp && base_len >= min_bases &&
          len * max_bases_per_snp >= base_len
        if (keep) {
          runs[[length(runs) + 1]] <- data.frame(
            first = idx[start + 1], last = idx[c], n_loci = len,
            n_het = n_het, n_miss = n_miss
          )
        }
        open <- FALSE
      }
      if (hit) {
        if (!open) {
          if (max_het == 0 && call_het) next
          open <- TRUE
          start <- c
          n_het <- 0
          n_miss <- 0
        }
        n_het <- n_het + call_het
        n_miss <- n_miss + is_miss[c + 1]
      }
    }
  }
  if (length(runs) == 0) {
    return(data.frame(
      first = integer(0), last = integer(0), n_loci = integer(0),
      n_het = integer(0), n_miss = integer(0)
    ))
  }
  do.call(rbind, runs)
}
