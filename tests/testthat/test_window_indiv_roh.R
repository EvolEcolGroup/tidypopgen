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
  # used to do)
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
      utils::getFromNamespace("slidingRuns", "detectRUNS")(
        genotypes[i, ],
        list(FID = groups[i], IID = x$id[i]),
        map,
        diff(map$bps),
        parameters
      )
    })
    runs <- dplyr::bind_rows(runs)
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

test_that("autoplot draws indiv_roh objects", {
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
  expect_s3_class(roh, "indiv_roh")
  # the attributes used by autoplot
  expect_identical(attr(roh, "individuals")$id, sheep_gt$id)
  expect_identical(
    as.character(attr(roh, "chromosome_length")$chromosome),
    c("2", "24")
  )
  for (type in c("runs", "stacked", "coverage", "froh", "violin",
                 "length_class")) {
    expect_s3_class(autoplot(roh, type = type), "ggplot")
  }
  expect_s3_class(
    autoplot(roh, type = "froh", froh_by = "chromosome"),
    "ggplot"
  )
  expect_s3_class(autoplot(roh, type = "runs", chromosome = "24"), "ggplot")
  expect_error(autoplot(roh, chromosome = "1"), "chromosomes with runs")
  # F_ROH is that of detectRUNS
  froh <- autoplot(roh, type = "froh")$data
  suppressMessages(
    ref <- detectRUNS::Froh_inbreeding(as_detectRUNS(roh), sheep_map)
  )
  expect_equal(froh$froh[match(ref$id, froh$id)], ref$Froh_genome)
})
