test_that("windows_indiv_roh works without backingfile update", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 1, 1, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0, 1),
    c(2, 2, 0, 0, 1, 1, 2)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:7),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2, 2)),
    position = as.integer(c(1000, 1500, 2000, 2500, 50000, 60000, 60500)),
    genetic_dist = as.double(rep(0, 7)),
    allele_ref = c("A", "T", "C", "G", "C", "T", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A", "A")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  suppressMessages(test_roh1 <- windows_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  # test_gt is not grouped, so there is no group column (but as_detectRUNS()
  # uses the id as the group)
  expect_false("group" %in% colnames(test_roh1))
  expect_equal(as_detectRUNS(test_roh1)$group, test_roh1$id)

  # remove one snp
  test_gt <- test_gt %>% select_loci(-7)

  # Check roh still works
  suppressMessages(test_roh2 <- windows_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  # the same runs (the length of chr2, in the attributes, has changed)
  expect_equal(test_roh1, test_roh2, ignore_attr = TRUE)

  # check that we get deprecation error if we use gt_roh_window
  suppressMessages(expect_warning(
    gt_roh_window(
      test_gt,
      window_size = 4,
      min_snp = 2,
      min_density = 1 / 500,
      max_gap = 5000,
      min_length_bps = 1000,
      threshold = 0.05
    ),
    "This is a soft-deprecated function,"
  ))
})

test_that("grouping method", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 1, 1, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0, 1),
    c(2, 2, 0, 0, 1, 1, 2)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:7),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2, 2)),
    position = as.integer(c(1000, 1500, 2000, 2500, 50000, 60000, 60500)),
    genetic_dist = as.double(rep(0, 7)),
    allele_ref = c("A", "T", "C", "G", "C", "T", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A", "A")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  test_gt <- test_gt %>% group_by(population)

  suppressMessages(test_roh1 <- windows_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  # test_gt is grouped, so group is grouping levels
  expect_identical(colnames(test_roh1)[1:2], c("group", "id"))
  expect_equal(test_roh1$group, dplyr::group_keys(test_gt) %>% pull(1))
  expect_equal(as_detectRUNS(test_roh1)$group, test_roh1$group)
})

test_that("results after subsetting", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d"),
    population = c("pop1", "pop1", "pop2", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 1, 1, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0, 1),
    c(2, 2, 2, 1, 0, 1, 2),
    c(2, 2, 0, 0, 1, 1, 2)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:7),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2, 2)),
    position = as.integer(c(1000, 1500, 2000, 2500, 50000, 60000, 60500)),
    genetic_dist = as.double(rep(0, 7)),
    allele_ref = c("A", "T", "C", "G", "C", "T", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A", "A")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  suppressMessages(test_roh_full <- windows_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  # subset individuals
  test_gt_sub_ind <- test_gt %>% filter(id %in% c("c", "d"))

  suppressMessages(test_roh_sub_ind <- windows_indiv_roh(
    test_gt_sub_ind,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  # results should be identical to the corresponding individuals in full results
  expect_true(all.equal(test_roh_full[c(2, 3), ],
    test_roh_sub_ind,
    check.attributes = FALSE
  ))
})

test_that("homozygotes for the alternate allele are homozygous", {
  test_genotypes <- rbind(
    rep(1, 6),
    rep(2, 6),
    rep(0, 6)
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
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
  # runs of heterozygosity are only found in the heterozygous individual
  suppressMessages(test_rohet <- windows_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    heterozygosity = TRUE,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))
  expect_equal(unique(test_rohet$id), "het")
  # runs of homozygosity are the same for both homozygotes
  suppressMessages(test_rohom <- windows_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    heterozygosity = FALSE,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))
  expect_equal(sort(unique(test_rohom$id)), c("hom_alt", "hom_ref"))
  run_cols <- c("chromosome", "n_loci", "start", "end", "length_bp")
  expect_equal(
    test_rohom[test_rohom$id == "hom_alt", run_cols],
    test_rohom[test_rohom$id == "hom_ref", run_cols],
    ignore_attr = TRUE
  )
})

test_that("windows_indiv_roh gives the same runs as detectRUNS", {
  skip_if_not_installed("detectRUNS")
  sheep_ped <- system.file("extdata", "Kijas2016_Sheep_subset.ped",
    package = "detectRUNS"
  )
  sheep_gt <- gen_tibble(sheep_ped, backingfile = tempfile(), quiet = TRUE)
  sheep_gt <- sheep_gt[c(1:6, 51:56), ] %>% group_by(population)
  # the runs of detectRUNS, one individual at a time (as windows_indiv_roh()
  # used to do), and one chromosome at a time (as windows restart at each
  # chromosome)
  ref_roh <- function(x, ...) {
    parameters <- list(
      windowSize = 15, threshold = 0.05, minSNP = 3, ROHet = FALSE,
      maxOppWindow = 1, maxMissWindow = 1, maxGap = 10^6,
      minLengthBps = 1000, minDensity = 1 / 1000, maxOppRun = NULL,
      maxMissRun = NULL
    )
    parameters[names(list(...))] <- list(...)
    map <- data.frame(
      Chrom = show_loci(x)$chromosome,
      SNP = show_loci(x)$name,
      bps = show_loci(x)$position
    )
    genotypes <- show_genotypes(x)
    genotypes[genotypes == 2] <- 0
    groups <- as.character(x$population)
    runs <- lapply(seq_len(nrow(x)), function(i) {
      lapply(unique(map$Chrom), function(chrom) {
        chrom_loci <- map$Chrom == chrom
        utils::getFromNamespace("slidingRuns", "detectRUNS")(
          genotypes[i, chrom_loci],
          list(FID = groups[i], IID = x$id[i]),
          map[chrom_loci, ],
          diff(map$bps[chrom_loci]),
          parameters
        )
      })
    })
    runs <- dplyr::bind_rows(unlist(runs, recursive = FALSE))
    row.names(runs) <- NULL
    runs
  }
  as_doubles <- function(runs) {
    for (col in c("nSNP", "from", "to", "lengthBps")) {
      runs[[col]] <- as.double(runs[[col]])
    }
    runs
  }
  suppressMessages(expected <- ref_roh(sheep_gt))
  expect_identical(
    as_doubles(as_detectRUNS(windows_indiv_roh(sheep_gt))),
    as_doubles(expected)
  )
  suppressMessages(
    expected <- ref_roh(sheep_gt,
      ROHet = TRUE, windowSize = 10, threshold = 0.1, maxOppWindow = 2,
      maxOppRun = 1, maxMissRun = 1
    )
  )
  expect_identical(
    as_doubles(as_detectRUNS(windows_indiv_roh(
      sheep_gt,
      heterozygosity = TRUE, window_size = 10, threshold = 0.1,
      max_opp_window = 2, max_opp_run = 1, max_miss_run = 1
    ))),
    as_doubles(expected)
  )
})

test_that("windows_indiv_roh checks the loci and uses several cores", {
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
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = data.frame(id = c("a", "b", "c")),
    quiet = TRUE,
    backingfile = tempfile()
  )
  roh_1 <- windows_indiv_roh(test_gt,
    window_size = 3, min_snp = 2, n_cores = 1
  )
  expect_true(nrow(roh_1) > 0)
  expect_identical(
    windows_indiv_roh(test_gt, window_size = 3, min_snp = 2, n_cores = 2),
    roh_1
  )
  # windows restart at each chromosome, so the runs of a chromosome are
  # those found on that chromosome alone (with windows that crossed into
  # chr2, which is homozygous for b, b would have a run on the last two loci
  # of chr1 with threshold = 0.5, and its run on chr1 would end with locus 6
  # with the default threshold)
  roh_half <- windows_indiv_roh(test_gt,
    window_size = 3, min_snp = 2, threshold = 0.5
  )
  roh_chr1 <- windows_indiv_roh(test_gt %>% select_loci(1:6),
    window_size = 3, min_snp = 2, threshold = 0.5
  )
  expect_equal(roh_half[roh_half$chromosome == "chr1", ], roh_chr1,
    ignore_attr = TRUE
  )
  expect_false("b" %in% roh_chr1$id)
  expect_identical(
    roh_1$last_locus[roh_1$id == "b" & roh_1$chromosome == "chr1"],
    5L
  )
  # with fewer loci than window_size, there are no runs
  roh_none <- windows_indiv_roh(test_gt, window_size = 13)
  expect_identical(nrow(roh_none), 0L)
  expect_identical(names(roh_none), names(roh_1))
  # the loci must be sorted by chromosome and position
  unsorted_gt <- test_gt %>% select_loci(c(1:3, 7:9, 4:6, 10:12))
  expect_error(
    windows_indiv_roh(unsorted_gt, window_size = 3),
    "sorted by chromosome and position"
  )
  unsorted_gt <- test_gt %>% select_loci(c(2, 1, 3:12))
  expect_error(
    windows_indiv_roh(unsorted_gt, window_size = 3),
    "sorted by chromosome and position"
  )
  # positions can not be missing
  loci <- show_loci(test_gt)
  loci$position[4] <- NA
  show_loci(test_gt) <- loci
  expect_error(
    windows_indiv_roh(test_gt, window_size = 3),
    "1 loci have a missing position"
  )
})

test_that("windows_indiv_roh returns tidy runs", {
  # individual a: a run on chr1 with one heterozygote and one missing
  # genotype; b: a run on chr2
  test_genotypes <- rbind(
    c(0, 0, 1, 0, NA, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1),
    c(1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 2, 2)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:16),
    chromosome = rep(c("chr1", "chr2"), each = 8),
    position = as.integer(rep(1:8 * 1000, 2)),
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = data.frame(id = c("a", "b"), pop = "p1", sex = c("F", "M")),
    quiet = TRUE,
    backingfile = tempfile()
  )
  roh <- windows_indiv_roh(test_gt,
    window_size = 4, min_snp = 3, max_opp_window = 1, max_miss_window = 1
  )
  expect_s3_class(roh, "tbl_df")
  expect_identical(names(roh), c(
    "id", "chromosome", "start", "end", "length_bp", "n_loci",
    "first_locus", "last_locus", "n_opp", "n_missing"
  ))
  expect_identical(levels(roh$chromosome), c("chr1", "chr2"))
  # the loci of each run are rows of the loci table
  loci <- show_loci(test_gt)
  expect_identical(roh$start, loci$position[roh$first_locus])
  expect_identical(roh$end, loci$position[roh$last_locus])
  expect_identical(roh$length_bp, roh$end - roh$start)
  expect_identical(
    as.character(roh$chromosome),
    as.character(loci$chromosome[roh$first_locus])
  )
  # the counts of opposite and missing genotypes are those of the run's loci
  geno <- show_genotypes(test_gt)
  for (r in seq_len(nrow(roh))) {
    run_loci <- seq(roh$first_locus[r], roh$last_locus[r])
    run_geno <- geno[match(roh$id[r], test_gt$id), run_loci]
    expect_identical(roh$n_loci[r], length(run_geno))
    expect_identical(roh$n_opp[r], sum(run_geno == 1, na.rm = TRUE))
    expect_identical(roh$n_missing[r], sum(is.na(run_geno)))
  }
  # with two grouping variables, one column for each
  roh_grouped <- windows_indiv_roh(test_gt %>% group_by(pop, sex),
    window_size = 4, min_snp = 3, max_opp_window = 1, max_miss_window = 1
  )
  expect_identical(names(roh_grouped)[1:3], c("pop", "sex", "id"))
  # the same runs (the individuals, in the attributes, now have groups)
  expect_equal(roh_grouped[-(1:2)], roh, ignore_attr = TRUE)
  expect_identical(
    as_detectRUNS(roh_grouped)$group,
    paste(roh_grouped$pop, roh_grouped$sex, sep = ".")
  )
  expect_error(as_detectRUNS(data.frame(id = 1)), "windows_indiv_roh")
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
  roh <- windows_indiv_roh(sim_gt %>% group_by(population))
  # each simulated run is found as a single run, which starts and ends at
  # most one window (15 loci) beyond it: loci outside the run that are
  # homozygous by chance are in homozygous windows that overlap the run
  for (i in seq_len(nrow(truth))) {
    overlaps <- roh$first_locus <= truth$last_locus[i] &
      roh$last_locus >= truth$first_locus[i]
    found <- roh[roh$id == truth$id[i] & overlaps, ]
    expect_identical(nrow(found), 1L)
    expect_identical(
      as.character(found$chromosome), as.character(truth$chromosome[i])
    )
    expect_true(found$first_locus <= truth$first_locus[i])
    expect_true(found$first_locus > truth$first_locus[i] - 15)
    expect_true(found$last_locus >= truth$last_locus[i])
    expect_true(found$last_locus < truth$last_locus[i] + 15)
  }
  # other runs are chance stretches of homozygous loci, much shorter than the
  # simulated runs (at least 60 loci)
  in_truth <- vapply(seq_len(nrow(roh)), function(r) {
    overlaps <- truth$first_locus <= roh$last_locus[r] &
      truth$last_locus >= roh$first_locus[r]
    any(truth$id == roh$id[r] & overlaps)
  }, logical(1))
  expect_true(all(roh$n_loci[!in_truth] < 30))
  # F_ROH (the proportion of the genome in runs) grows with inbreeding
  froh <- autoplot(roh, type = "froh")$data
  mean_froh <- tapply(froh$froh, froh$group, mean)
  expect_true(mean_froh[["outbred"]] < 0.01)
  expect_true(mean_froh[["outbred"]] < mean_froh[["moderate"]])
  expect_true(mean_froh[["moderate"]] < mean_froh[["inbred"]])
})

test_that("windows_indiv_roh deals with the edge cases of the simulated runs", {
  sim_gt <- roh_sim_gt()
  loci <- show_loci(sim_gt)
  roh <- windows_indiv_roh(sim_gt)
  # inbred_01: runs at the end of chromosome 1 and at the start of chromosome
  # 2 are two runs, which end and start at the last and first locus of their
  # chromosome
  runs_01 <- roh[roh$id == "inbred_01" & roh$n_loci > 30, ]
  expect_identical(as.character(runs_01$chromosome), c("1", "2"))
  expect_identical(runs_01$last_locus[1], max(which(loci$chromosome == "1")))
  expect_identical(runs_01$first_locus[2], min(which(loci$chromosome == "2")))
  # inbred_02: the run on chromosome 3 holds a heterozygous genotype (locus
  # 2450) and two missing genotypes (loci 2400 and 2500), but no window holds
  # more than one of each, so it is a single run
  has_2450 <- roh$first_locus <= 2450 & roh$last_locus >= 2450
  run_02 <- roh[roh$id == "inbred_02" & has_2450, ]
  expect_identical(nrow(run_02), 1L)
  expect_true(run_02$n_opp >= 1)
  expect_identical(run_02$n_missing, 2L)
  # it is dropped if runs can not have heterozygous genotypes
  roh_no_opp <- windows_indiv_roh(sim_gt[sim_gt$id == "inbred_02", ],
    max_opp_run = 0
  )
  expect_false(any(
    roh_no_opp$first_locus <= 2450 & roh_no_opp$last_locus >= 2450
  ))
  # and it is split around them if windows can not hold any
  roh_strict <- windows_indiv_roh(sim_gt[sim_gt$id == "inbred_02", ],
    max_opp_window = 0, max_miss_window = 0
  )
  roh_strict <- roh_strict[roh_strict$n_loci > 30, ]
  expect_identical(roh_strict$last_locus[1:3], c(2399L, 2449L, 2499L))
  expect_identical(roh_strict$first_locus[2:4], c(2401L, 2451L, 2501L))
  # inbred_03: a gap larger than max_gap makes the windows that span it
  # non-homozygous, but, as in detectRUNS, it does not split the run, as the
  # loci on each side of the gap are in homozygous windows
  gap <- which(diff(loci$position) > 10^6)
  spans_gap <- roh$first_locus <= gap & roh$last_locus > gap
  run_03 <- roh[roh$id == "inbred_03" & spans_gap, ]
  expect_identical(nrow(run_03), 1L)
})

test_that("autoplot draws indiv_roh objects", {
  sim_gt <- roh_sim_gt() %>% group_by(population)
  roh <- windows_indiv_roh(sim_gt)
  expect_s3_class(roh, "indiv_roh")
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
  sheep_ped <- system.file("extdata", "Kijas2016_Sheep_subset.ped",
    package = "detectRUNS"
  )
  sheep_map <- system.file("extdata", "Kijas2016_Sheep_subset.map",
    package = "detectRUNS"
  )
  sheep_gt <- gen_tibble(sheep_ped, backingfile = tempfile(), quiet = TRUE)
  sheep_gt <- sheep_gt[c(1:6, 51:56), ] %>% group_by(population)
  roh <- windows_indiv_roh(sheep_gt)
  froh <- autoplot(roh, type = "froh")$data
  suppressMessages(
    ref <- detectRUNS::Froh_inbreeding(as_detectRUNS(roh), sheep_map)
  )
  expect_equal(froh$froh[match(ref$id, froh$id)], ref$Froh_genome)
})
