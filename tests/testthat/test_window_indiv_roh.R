# small parameters, for the small examples below
small_par <- list(
  window_snp = 3, min_snp = 2, min_kb = 0, max_kb_per_snp = 50,
  max_gap_kb = 5
)
small_roh <- function(x, ...) {
  args <- utils::modifyList(small_par, list(...))
  do.call(windows_indiv_roh, c(list(x), args))
}

# a small gen_tibble: a and b have runs on chr1, c on chr2
small_gt <- function() {
  test_genotypes <- rbind(
    c(0, 0, 0, 0, 1, 0, 0, 0, NA, 0, 2, 2),
    c(0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0),
    c(2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:12),
    chromosome = rep(c("chr1", "chr2"), each = 6),
    position = as.integer(rep(1:6 * 1000, 2)),
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
  gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = data.frame(
      id = c("a", "b", "c"), population = c("pop1", "pop1", "pop2"),
      sex = c("F", "M", "F")
    ),
    quiet = TRUE,
    backingfile = tempfile()
  )
}

test_that("windows_indiv_roh returns tidy runs", {
  test_gt <- small_gt()
  roh <- small_roh(test_gt, n_cores = 1)
  expect_s3_class(roh, "indiv_roh")
  expect_s3_class(roh, "tbl_df")
  expect_identical(names(roh), c(
    "id", "chromosome", "start", "end", "length_bp", "n_loci",
    "first_locus", "last_locus", "n_het", "n_missing"
  ))
  expect_identical(levels(roh$chromosome), c("chr1", "chr2"))
  expect_true(nrow(roh) > 0)
  # the loci of each run are rows of the loci table
  loci <- show_loci(test_gt)
  expect_identical(roh$start, loci$position[roh$first_locus])
  expect_identical(roh$end, loci$position[roh$last_locus])
  # lengths are those of PLINK 1.9, or of PLINK 1.07 (one less)
  expect_identical(roh$length_bp, roh$end - roh$start + 1L)
  roh_old <- small_roh(test_gt, subtract_1_from_lengths = TRUE)
  expect_identical(roh_old$length_bp, roh_old$end - roh_old$start)
  expect_identical(
    as.character(roh$chromosome),
    as.character(loci$chromosome[roh$first_locus])
  )
  # the counts of heterozygous and missing genotypes are those of the run
  geno <- show_genotypes(test_gt)
  for (r in seq_len(nrow(roh))) {
    run_loci <- seq(roh$first_locus[r], roh$last_locus[r])
    run_geno <- geno[match(roh$id[r], test_gt$id), run_loci]
    expect_identical(roh$n_loci[r], length(run_geno))
    expect_identical(roh$n_het[r], sum(run_geno == 1, na.rm = TRUE))
    expect_identical(roh$n_missing[r], sum(is.na(run_geno)))
  }
  # the same runs with several cores
  expect_identical(small_roh(test_gt, n_cores = 2), roh)
  # an ungrouped gen_tibble has no group column (as_detectRUNS() uses the
  # id as the group)
  expect_equal(as_detectRUNS(roh)$group, roh$id)
  # a single grouping variable gives a group column
  roh_pop <- small_roh(test_gt %>% group_by(population))
  expect_identical(names(roh_pop)[1:2], c("group", "id"))
  expect_identical(
    roh_pop$group,
    test_gt$population[match(roh_pop$id, test_gt$id)]
  )
  # with two grouping variables, one column for each
  roh_grouped <- small_roh(test_gt %>% group_by(population, sex))
  expect_identical(names(roh_grouped)[1:3], c("population", "sex", "id"))
  # the same runs (the individuals, in the attributes, now have groups)
  expect_equal(roh_grouped[-(1:2)], roh, ignore_attr = TRUE)
  expect_identical(
    as_detectRUNS(roh_grouped)$group,
    paste(roh_grouped$population, roh_grouped$sex, sep = ".")
  )
  expect_error(as_detectRUNS(data.frame(id = 1)), "windows_indiv_roh")
  # the runs of a subset of individuals are theirs in the full dataset
  roh_sub <- small_roh(test_gt %>% filter(id %in% c("b", "c")))
  expect_equal(roh_sub, roh[roh$id %in% c("b", "c"), ], ignore_attr = TRUE)
})

test_that("windows_indiv_roh checks its arguments and the loci", {
  test_gt <- small_gt()
  # windows restart at each chromosome: with fewer loci than window_snp,
  # there are no runs
  roh_none <- small_roh(test_gt, window_snp = 7)
  expect_identical(nrow(roh_none), 0L)
  expect_identical(names(roh_none), names(small_roh(test_gt)))
  # the arguments of the version based on detectRUNS give an error
  expect_error(
    windows_indiv_roh(test_gt, window_size = 3, threshold = 0.05),
    "window_size, threshold are no longer used"
  )
  expect_error(windows_indiv_roh(test_gt, foo = 1), "must be empty")
  expect_error(small_roh(test_gt, window_snp = 0), "window_snp")
  expect_error(small_roh(test_gt, max_het = -1), "max_het")
  expect_error(small_roh(test_gt, window_threshold = 2), "window_threshold")
  # the loci must be sorted by chromosome and position
  unsorted_gt <- test_gt %>% select_loci(c(1:3, 7:9, 4:6, 10:12))
  expect_error(small_roh(unsorted_gt), "sorted by chromosome and position")
  unsorted_gt <- test_gt %>% select_loci(c(2, 1, 3:12))
  expect_error(small_roh(unsorted_gt), "sorted by chromosome and position")
  # positions can not be missing
  loci <- show_loci(test_gt)
  loci$position[4] <- NA
  show_loci(test_gt) <- loci
  expect_error(small_roh(test_gt), "1 loci have a missing position")
})

test_that("homozygotes for the alternate allele are homozygous", {
  test_gt <- gen_tibble(
    x = rbind(rep(1, 6), rep(2, 6), rep(0, 6)),
    loci = data.frame(
      name = paste0("rs", 1:6),
      chromosome = "chr1",
      position = as.integer(1:6 * 500 + 500),
      genetic_dist = as.double(rep(0, 6)),
      allele_ref = "A",
      allele_alt = "T"
    ),
    indiv_meta = data.frame(
      id = c("het", "hom_alt", "hom_ref"),
      population = "pop1"
    ),
    quiet = TRUE,
    backingfile = tempfile()
  )
  roh <- small_roh(test_gt)
  expect_identical(roh$id, c("hom_alt", "hom_ref"))
  run_cols <- c("chromosome", "n_loci", "start", "end", "length_bp")
  expect_equal(roh[1, run_cols], roh[2, run_cols], ignore_attr = TRUE)
})

test_that("windows_indiv_roh follows the algorithm of PLINK --homozyg", {
  # compare with an R transcription of PLINK's code (helper_plink_homozyg.R)
  # on random genotypes, with parameters that make its edge cases common
  # (short chromosomes, gaps, heterozygous genotypes in runs)
  set.seed(42)
  for (it in 1:30) {
    n_chrom <- sample(1:3, 1)
    n_per <- sample(1:30, n_chrom, replace = TRUE)
    m <- sum(n_per)
    chrom <- rep(seq_len(n_chrom), n_per)
    pos <- unlist(lapply(n_per, function(k) {
      cumsum(sample(c(1:2000, 5000:9000), k, replace = TRUE))
    }))
    geno <- matrix(
      sample(c(0, 1, 2, NA), 3 * m,
        replace = TRUE,
        prob = c(0.4, 0.15, 0.4, 0.05)
      ),
      nrow = 3
    )
    par <- list(
      window_snp = sample(1:6, 1), window_het = sample(0:1, 1),
      window_missing = sample(0:1, 1),
      window_threshold = sample(c(0, 0.05, 0.2, 0.5, 1), 1),
      min_snp = sample(1:4, 1), min_kb = sample(c(0, 3), 1),
      max_kb_per_snp = sample(c(2, 50), 1), max_gap_kb = sample(c(3, 1000), 1),
      max_het = sample(c(0, 1, Inf), 1),
      subtract_1_from_lengths = sample(c(TRUE, FALSE), 1)
    )
    test_gt <- gen_tibble(geno,
      indiv_meta = data.frame(id = c("a", "b", "c"), population = "p"),
      loci = data.frame(
        name = paste0("s", seq_len(m)), chromosome = chrom,
        position = as.integer(pos), genetic_dist = 0, allele_ref = "A",
        allele_alt = "T"
      ),
      quiet = TRUE,
      backingfile = tempfile()
    )
    roh <- do.call(windows_indiv_roh, c(list(test_gt, n_cores = 1), par))
    expected <- do.call(rbind, lapply(1:3, function(i) {
      runs <- do.call(plink_homozyg_r, c(list(geno[i, ], chrom, pos), par))
      runs$id <- rep(test_gt$id[i], nrow(runs))
      runs
    }))
    expect_identical(roh$id, expected$id)
    expect_identical(roh$first_locus, as.integer(expected$first))
    expect_identical(roh$last_locus, as.integer(expected$last))
    expect_identical(roh$n_loci, as.integer(expected$n_loci))
    expect_identical(roh$n_het, as.integer(expected$n_het))
    expect_identical(roh$n_missing, as.integer(expected$n_miss))
  }
})

# the simulated dataset of inst/extdata/roh, with known runs of homozygosity
# (see data-raw/create_roh_dataset.R)
roh_sim_gt <- function() {
  gen_tibble(
    system.file("extdata", "roh", "roh_sim.bed", package = "tidypopgen"),
    backingfile = tempfile(),
    quiet = TRUE
  )
}
roh_sim_truth <- function() {
  utils::read.csv(
    system.file("extdata", "roh", "roh_sim_truth.csv", package = "tidypopgen")
  )
}

test_that("windows_indiv_roh finds the simulated runs", {
  sim_gt <- roh_sim_gt()
  truth <- roh_sim_truth()
  # the default parameters, those of PLINK
  roh <- windows_indiv_roh(sim_gt %>% group_by(population))
  # each simulated run is found (the run of inbred_03 across a large gap is
  # found as two runs, see below), and its ends are within a few loci of the
  # true ones: windows are hits only with at most one heterozygous genotype,
  # which a window that overlaps the end of a run rarely has, so runs can be
  # a few loci shorter; and loci outside the run that are homozygous by
  # chance can be added
  for (i in seq_len(nrow(truth))) {
    overlaps <- roh$first_locus <= truth$last_locus[i] &
      roh$last_locus >= truth$first_locus[i]
    found <- roh[roh$id == truth$id[i] & overlaps, ]
    expect_identical(nrow(found), if (truth$id[i] == "inbred_03") 2L else 1L)
    expect_identical(
      as.character(found$chromosome[1]), as.character(truth$chromosome[i])
    )
    expect_true(abs(found$first_locus[1] - truth$first_locus[i]) < 10)
    expect_true(abs(found$last_locus[nrow(found)] - truth$last_locus[i]) < 10)
  }
  # there are no other runs (with windows of 50 loci, unrelated haplotypes
  # never look like a run)
  expect_identical(nrow(roh), nrow(truth) + 1L)
  # F_ROH (the proportion of the genome in runs) grows with inbreeding
  froh <- autoplot(roh, type = "froh")$data
  mean_froh <- tapply(froh$froh, froh$group, mean)
  expect_identical(mean_froh[["outbred"]], 0)
  expect_true(mean_froh[["moderate"]] < mean_froh[["inbred"]])
})

test_that("windows_indiv_roh deals with the edge cases of the simulated runs", {
  sim_gt <- roh_sim_gt()
  loci <- show_loci(sim_gt)
  roh <- windows_indiv_roh(sim_gt)
  # inbred_01: runs at the end of chromosome 1 and at the start of chromosome
  # 2 are two runs, which end and start at the last and first locus of their
  # chromosome (or very near it)
  runs_01 <- roh[roh$id == "inbred_01", ]
  expect_identical(as.character(runs_01$chromosome), c("1", "2"))
  expect_identical(runs_01$last_locus[1], max(which(loci$chromosome == "1")))
  expect_true(runs_01$first_locus[2] - min(which(loci$chromosome == "2")) < 3)
  # inbred_02: the run on chromosome 3 holds a heterozygous genotype (locus
  # 2450) and missing genotypes (loci 2400 and 2500), but no window holds
  # more than one heterozygous genotype, so it is a single run
  inbred_02 <- sim_gt[sim_gt$id == "inbred_02", ]
  run_02 <- roh[roh$id == "inbred_02" & roh$chromosome == "3", ]
  expect_identical(nrow(run_02), 1L)
  expect_identical(run_02$n_het, 1L)
  expect_true(run_02$n_missing >= 2)
  # with no heterozygous genotypes allowed in a run, it is split at the
  # heterozygous genotype (which is in neither run); each part is shorter
  # than 1 Mb
  roh_no_het <- windows_indiv_roh(inbred_02, max_het = 0, min_kb = 500)
  roh_no_het <- roh_no_het[roh_no_het$chromosome == "3", ]
  expect_identical(roh_no_het$last_locus[1], 2449L)
  expect_identical(roh_no_het$first_locus[2], 2451L)
  expect_identical(roh_no_het$n_het, c(0L, 0L))
  # with windows that can not hold it, it is split too (the windows that
  # include it are not hits, so the loci around it are not in a run)
  roh_strict <- windows_indiv_roh(inbred_02, window_het = 0, min_kb = 500)
  roh_strict <- roh_strict[roh_strict$chromosome == "3", ]
  expect_identical(nrow(roh_strict), 2L)
  expect_true(roh_strict$last_locus[1] < 2450)
  expect_true(roh_strict$first_locus[2] > 2450)
  # inbred_03: a gap larger than max_gap_kb splits its run on chromosome 2
  # (windows span the gap, so each part extends to the gap)
  gap <- which(diff(loci$position) > 10^6)
  runs_03 <- roh[roh$id == "inbred_03", ]
  expect_identical(nrow(runs_03), 2L)
  expect_identical(runs_03$last_locus[1], gap)
  expect_identical(runs_03$first_locus[2], gap + 1L)
  # unless the largest gap is longer
  runs_03 <- windows_indiv_roh(sim_gt[sim_gt$id == "inbred_03", ],
    max_gap_kb = 2000
  )
  expect_identical(nrow(runs_03), 1L)
})

test_that("autoplot draws indiv_roh objects", {
  sim_gt <- roh_sim_gt() %>% group_by(population)
  roh <- windows_indiv_roh(sim_gt)
  # the attributes used by autoplot
  expect_identical(attr(roh, "individuals")$id, sim_gt$id)
  expect_identical(
    as.character(attr(roh, "chromosome_length")$chromosome),
    c("1", "2", "3")
  )
  for (type in c("runs", "stacked", "coverage", "froh", "violin",
                 "length_class")) {
    expect_s3_class(autoplot(roh, type = type), "ggplot")
  }
  expect_s3_class(
    autoplot(roh, type = "froh", froh_by = "chromosome"),
    "ggplot"
  )
  expect_s3_class(autoplot(roh, type = "runs", chromosome = "3"), "ggplot")
  expect_error(autoplot(roh, chromosome = "4"), "chromosomes with runs")
})

test_that("F_ROH of autoplot is that of detectRUNS", {
  skip_if_not_installed("detectRUNS")
  sim_gt <- roh_sim_gt() %>% group_by(population)
  roh <- windows_indiv_roh(sim_gt)
  # the map of detectRUNS
  map_file <- tempfile(fileext = ".map")
  loci <- show_loci(sim_gt)
  utils::write.table(
    data.frame(loci$chromosome, loci$name, 0, loci$position),
    map_file,
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  froh <- autoplot(roh, type = "froh")$data
  suppressMessages(
    ref <- detectRUNS::Froh_inbreeding(as_detectRUNS(roh), map_file)
  )
  expect_equal(froh$froh[match(ref$id, froh$id)], ref$Froh_genome)
})

test_that("gt_roh_window is the deprecated version based on detectRUNS", {
  skip_if_not_installed("detectRUNS")
  sim_gt <- roh_sim_gt()[c(1, 11, 21:23), ] %>% group_by(population)
  expect_warning(
    suppressMessages(roh_old <- gt_roh_window(sim_gt)),
    "gt_roh_window\\(\\) is deprecated"
  )
  # the format of detectRUNS
  expect_identical(
    names(roh_old),
    c("group", "id", "chrom", "nSNP", "from", "to", "lengthBps")
  )
  expect_true(all(roh_old$group %in% c("outbred", "moderate", "inbred")))
  # and its runs of heterozygosity
  expect_warning(suppressMessages(
    rohet <- gt_roh_window(sim_gt, heterozygosity = TRUE, window_size = 10)
  ))
  expect_identical(names(rohet), names(roh_old))
})
