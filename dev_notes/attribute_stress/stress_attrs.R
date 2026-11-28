# Stress test: do "ancestry_labels" (extended phased codes, #245) and
# "chr_types" (declared chromosome types, sex_chromosomes) survive any order
# of operations? Run from the package root of the sex_chromosomes worktree.
pkgload::load_all(quiet = TRUE)
suppressPackageStartupMessages(library(dplyr))
set.seed(2026)

n_ind <- 12
n_loci <- 40
labels <- c("anc_a", "anc_b")
chr_names <- rep(c("1", "2", "23", "26"), each = n_loci / 4)
# PLINK's human codes: 23 is X, 26 is MT
declared_types <- gt_set_chromosome_types # (alias, for readability)

extended_code <- function(k) {
  n_blocks <- 1 + k^2
  c(rep(c(0, 1, 2, NA, 1), n_blocks), rep(NA, 256 - 5 * n_blocks))
}

make_haps <- function(n_ind, n_loci) {
  freq <- stats::runif(n_loci, 0.05, 0.95)
  haps <- matrix(
    stats::rbinom(2 * n_ind * n_loci, 1, rep(freq, each = 2 * n_ind)),
    nrow = 2 * n_ind
  )
  missing <- matrix(stats::runif(n_ind * n_loci) < 0.05, nrow = n_ind)
  haps[missing[rep(seq_len(n_ind), each = 2), , drop = FALSE]] <- NA
  storage.mode(haps) <- "integer"
  haps
}

make_ext <- function(haps, prefix, labels) {
  n <- nrow(haps) / 2
  indiv <- data.frame(
    id = paste0(prefix, seq_len(n)),
    population = rep(c("p1", "p2", "p3"), length.out = n)
  )
  loci <- data.frame(
    name = paste0("rs", seq_len(ncol(haps))),
    chromosome = chr_names,
    position = as.integer(seq_len(ncol(haps)) * 100),
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
  x <- gen_tibble(haps,
    phased = TRUE, indiv_meta = indiv, loci = loci,
    backingfile = tempfile("ext_"), quiet = TRUE
  )
  fbm <- .gt_get_fbm(x)
  bytes <- fbm$copy(code = as.numeric(0:255))[, ]
  fbm[] <- as.raw(bytes + 5 * sample(0:4, length(bytes), replace = TRUE))
  fbm$code256 <- extended_code(2)
  fbm$save()
  attr(x$genotypes, "fbm_md5sum") <- tools::md5sum(fbm$rds)
  attr(x$genotypes, "ancestry_labels") <- labels
  x
}

haps_a <- make_haps(n_ind, n_loci)
haps_b <- make_haps(n_ind, n_loci)
x0 <- make_ext(haps_a, "a", labels)
other <- make_ext(haps_b, "b", labels)
# the reference haplotypes, by id and locus name
ref_haps <- rbind(haps_a, haps_b)
rownames(ref_haps) <- paste0(
  rep(c(paste0("a", 1:n_ind), paste0("b", 1:n_ind)), each = 2),
  c("_1", "_2")
)
colnames(ref_haps) <- paste0("rs", seq_len(n_loci))

expected_types <- attr(
  gt_set_chromosome_types(x0, chr_set = "human")$genotypes, "chr_types"
)

# the operations; each gets the current gen_tibble and returns a new one
ops <- list(
  declare = function(x) gt_set_chromosome_types(x, chr_set = "human"),
  filter_ind = function(x) {
    keep <- sample(x$id, max(4, ceiling(0.8 * nrow(x))))
    filter(x, id %in% keep)
  },
  bracket_rows = function(x) x[sample(nrow(x)), ],
  arrange = function(x) arrange(x, desc(id)),
  group_mutate = function(x) {
    x %>%
      group_by(population) %>%
      mutate(n_pop = n()) %>%
      ungroup()
  },
  select_loci = function(x) {
    nm <- loci_names(x)
    # drop a few loci, but keep some of each chromosome, and shuffle them
    keep <- unlist(lapply(split(nm, show_loci(x)$chromosome), function(v) {
      sample(v, max(2, length(v) - 2))
    }))
    select_loci(x, all_of(sample(keep)))
  },
  select_loci_if = function(x) {
    select_loci_if(x, stats::runif(count_loci(x)) > 0.1 |
      loci_chromosomes(genotypes) %in% c("23", "26"))
  },
  order_loci = function(x) gt_order_loci(x, quiet = TRUE),
  update_bk = function(x) {
    gt_update_backingfile(x, backingfile = tempfile("upd_"), quiet = TRUE)
  },
  save_load = function(x) {
    f <- tempfile("saved_")
    gt_save(x, file_name = f, quiet = TRUE)
    gt_load(f)
  },
  rbind = function(x) {
    # the other dataset, on the same loci, with the same declaration (if any)
    o <- select_loci(other, all_of(loci_names(x)))
    o <- filter(o, !id %in% x$id)
    if (!is.null(attr(x$genotypes, "chr_types"))) {
      o <- gt_set_chromosome_types(o, chr_set = "human")
    }
    suppressMessages(rbind(x, o,
      quiet = TRUE,
      backingfile = tempfile("rb_")
    ))
  },
  rbind_undeclared = function(x) {
    # the other dataset declares nothing: x's declaration should be kept
    o <- select_loci(other, all_of(loci_names(x)))
    o <- filter(o, !id %in% x$id)
    suppressMessages(rbind(x, o,
      quiet = TRUE,
      backingfile = tempfile("rbu_")
    ))
  },
  analysis = function(x) {
    # autosome-only analyses must not change x
    invisible(suppressMessages(indiv_het_obs(x)))
    invisible(suppressMessages(loci_maf(x)))
    x
  }
)

check <- function(x, declared, label) {
  problems <- character(0)
  if (!is_phased(x)) problems <- c(problems, "not phased")
  if (!identical(attr(x$genotypes, "ancestry_labels", exact = TRUE), labels)) {
    problems <- c(problems, "ancestry_labels lost or changed")
  }
  if (!identical(.gt_get_fbm(x)$code256, extended_code(2))) {
    problems <- c(problems, "code256 is no longer the extended code")
  }
  ct <- attr(x$genotypes, "chr_types")
  if (declared && !identical(ct, expected_types)) {
    problems <- c(problems, "chr_types lost or changed")
  }
  if (!declared && !is.null(ct)) {
    problems <- c(problems, "chr_types appeared")
  }
  # the types of the loci
  lt <- as.character(loci_chromosome_types(x$genotypes))
  chr <- as.character(show_loci(x)$chromosome)
  want <- if (declared) {
    c("1" = "autosome", "2" = "autosome", "23" = "X", "26" = "MT")[chr]
  } else {
    rep("autosome", length(chr))
  }
  if (!identical(unname(lt), unname(want))) {
    problems <- c(problems, "loci_chromosome_types wrong")
  }
  # the genotypes and the phase
  h <- show_haplotypes(x)
  rows <- paste0(rep(as.character(x$id), each = 2), c("_1", "_2"))
  want_h <- ref_haps[rows, loci_names(x), drop = FALSE]
  if (!isTRUE(all.equal(unname(h), unname(want_h)))) {
    problems <- c(problems, "haplotypes differ from the reference")
  }
  g <- show_genotypes(x)
  want_g <- want_h[c(TRUE, FALSE), , drop = FALSE] +
    want_h[c(FALSE, TRUE), , drop = FALSE]
  if (!isTRUE(all.equal(unname(g), unname(want_g)))) {
    problems <- c(problems, "genotypes differ from the reference")
  }
  if (length(problems)) paste0(label, ": ", paste(problems, collapse = "; "))
}

run_sequence <- function(seq_ops) {
  x <- x0
  declared <- FALSE
  for (i in seq_along(seq_ops)) {
    op <- seq_ops[i]
    x <- tryCatch(ops[[op]](x), error = function(e) e)
    if (inherits(x, "error")) {
      return(paste0(
        paste(seq_ops[1:i], collapse = " > "), ": ERROR ",
        conditionMessage(x)
      ))
    }
    if (op == "declare") declared <- TRUE
    res <- check(x, declared, paste(seq_ops[1:i], collapse = " > "))
    if (!is.null(res)) {
      return(res)
    }
  }
  NULL
}

op_names <- names(ops)
# every ordered pair, then random sequences of 6 operations
pairs <- expand.grid(a = op_names, b = op_names, stringsAsFactors = FALSE)
seqs <- c(
  lapply(seq_len(nrow(pairs)), function(i) unlist(pairs[i, ])),
  lapply(1:300, function(i) sample(op_names, 6, replace = TRUE))
)
# rbind twice can run out of new ids: only one rbind per sequence
seqs <- lapply(seqs, function(s) {
  r <- which(s %in% c("rbind", "rbind_undeclared"))
  if (length(r) > 1) s <- s[-r[-1]]
  s
})
results <- lapply(seqs, run_sequence)
failures <- unique(unlist(results))
cat(length(seqs), "sequences,", sum(lengths(seqs)), "operations;",
  sum(!vapply(results, is.null, logical(1))), "failed\n")
if (length(failures)) {
  # group by the failing step, to see the pattern
  cat(head(failures, 40), sep = "\n")
}
