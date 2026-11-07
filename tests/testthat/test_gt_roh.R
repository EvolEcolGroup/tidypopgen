# the runs of the simulated dataset (see data-raw/create_roh_dataset.R)
sim_gt <- gen_tibble(
  system.file("extdata", "roh", "roh_sim.bed", package = "tidypopgen"),
  backingfile = tempfile(),
  quiet = TRUE
) %>% group_by(population)
sim_roh <- windows_indiv_roh(sim_gt)

test_that("gt_roh_froh computes F_ROH", {
  froh <- gt_roh_froh(sim_roh)
  expect_identical(names(froh), c("group", "id", "length_bp", "froh"))
  # every individual, in the order of the gen_tibble, with 0 if no runs
  expect_identical(froh$id, sim_gt$id)
  expect_true(all(froh$froh[froh$group == "outbred"] == 0))
  # the total length of the runs over the length of the genome (the sum of
  # the largest position of each chromosome)
  loci <- show_loci(sim_gt)
  genome <- sum(tapply(loci$position, loci$chromosome, max))
  total <- tapply(sim_roh$length_bp, sim_roh$id, sum)
  expect_equal(
    froh$froh[match(names(total), froh$id)],
    as.numeric(total) / genome
  )
  # the same F_ROH as autoplot()
  plotted <- autoplot(sim_roh, type = "froh")$data
  expect_equal(plotted$froh[match(froh$id, plotted$id)], froh$froh)
  # a given length of the genome
  expect_equal(
    gt_roh_froh(sim_roh, genome_length = 3e9)$froh,
    froh$length_bp / 3e9
  )
  # from the runs of at least some lengths: one row per individual and
  # length, and F_ROH decreases with the length
  froh_class <- gt_roh_froh(sim_roh, min_length_mb = c(0, 1, 2, 4))
  expect_identical(names(froh_class)[1:3], c("group", "id", "min_length_mb"))
  expect_identical(froh_class$id, rep(sim_gt$id, each = 4))
  expect_equal(froh_class$froh[froh_class$min_length_mb == 0], froh$froh)
  expect_true(all(tapply(froh_class$froh, froh_class$id, function(x) {
    !is.unsorted(rev(x))
  })))
  long <- sim_roh[sim_roh$length_bp >= 2e6, ]
  expect_equal(
    sum(froh_class$length_bp[froh_class$min_length_mb == 2]),
    sum(as.numeric(long$length_bp))
  )
  # on each chromosome: the totals add up to those of the genome
  froh_chrom <- gt_roh_froh(sim_roh, by = "chromosome")
  expect_identical(nrow(froh_chrom), nrow(sim_gt) * 3L)
  expect_equal(
    as.numeric(tapply(froh_chrom$length_bp, froh_chrom$id, sum)[froh$id]),
    froh$length_bp
  )
  # errors
  expect_error(gt_roh_froh(tibble::as_tibble(sim_roh)), "windows_indiv_roh")
  expect_error(gt_roh_froh(sim_roh, min_length_mb = "a"), "min_length_mb")
  expect_error(gt_roh_froh(sim_roh, genome_length = -1), "genome_length")
})

test_that("gt_roh_froh gives the F_ROH of detectRUNS", {
  skip_if_not_installed("detectRUNS")
  map_file <- tempfile(fileext = ".map")
  loci <- show_loci(sim_gt)
  utils::write.table(
    data.frame(loci$chromosome, loci$name, 0, loci$position),
    map_file,
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  runs <- as_detectRUNS(sim_roh)
  froh <- gt_roh_froh(sim_roh)
  suppressMessages(ref <- detectRUNS::Froh_inbreeding(runs, map_file))
  expect_equal(froh$froh[match(ref$id, froh$id)], ref$Froh_genome)
  # by length class (with Class = 1, detectRUNS uses runs of at least 0, 1,
  # 2, 4 and 8 Mb, with no column for lengths without runs)
  froh_class <- gt_roh_froh(sim_roh, min_length_mb = c(1, 2))
  suppressMessages(utils::capture.output(
    ref_class <- detectRUNS::Froh_inbreedingClass(runs, map_file, Class = 1)
  ))
  for (min_length in c(1, 2)) {
    ref_col <- ref_class[[paste0("Froh_Class_", min_length)]]
    ref_col[is.na(ref_col)] <- 0
    ours <- froh_class[froh_class$min_length_mb == min_length, ]
    expect_equal(ours$froh[match(ref_class$id, ours$id)], ref_col)
  }
})

test_that("gt_roh_loci counts the individuals in runs at each locus", {
  loci <- show_loci(sim_gt)
  roh_loci <- gt_roh_loci(sim_roh, loci)
  expect_identical(
    names(roh_loci),
    c("group", "name", "chromosome", "position", "n_roh", "n_indiv",
      "prop_roh")
  )
  expect_identical(nrow(roh_loci), 3L * nrow(loci))
  # brute force: the runs of each group that include each locus
  for (g in c("inbred", "outbred")) {
    group_runs <- sim_roh[sim_roh$group == g, ]
    expected <- vapply(seq_len(nrow(loci)), function(k) {
      sum(group_runs$first_locus <= k & group_runs$last_locus >= k)
    }, integer(1))
    expect_identical(roh_loci$n_roh[roh_loci$group == g], expected)
  }
  expect_identical(unique(roh_loci$n_indiv), 10L)
  expect_equal(roh_loci$prop_roh, roh_loci$n_roh / 10)
  # over all individuals
  all_loci <- gt_roh_loci(sim_roh, loci, by_group = FALSE)
  expect_identical(names(all_loci)[1], "name")
  expect_identical(
    all_loci$n_roh,
    as.integer(rowSums(matrix(roh_loci$n_roh, ncol = 3)))
  )
  expect_identical(unique(all_loci$n_indiv), 30L)
  # the loci table must be that of the runs
  expect_error(gt_roh_loci(sim_roh, loci[-1, ]), "used to find the runs")
  expect_error(gt_roh_loci(sim_roh, data.frame(a = 1)), "loci table")
})

test_that("gt_roh_islands finds the regions in runs in many individuals", {
  loci <- show_loci(sim_gt)
  roh_loci <- gt_roh_loci(sim_roh, loci)
  islands <- gt_roh_islands(sim_roh, loci, threshold = 0.2)
  expect_identical(names(islands), c(
    "group", "chromosome", "start", "end", "length_bp", "n_loci",
    "first_locus", "last_locus", "max_prop", "mean_prop"
  ))
  expect_true(nrow(islands) > 0)
  # the loci of the islands are those with a proportion of at least 0.2, and
  # islands are separated by other loci or by the end of a chromosome
  for (g in unique(islands$group)) {
    g_islands <- islands[islands$group == g, ]
    in_island <- unlist(Map(seq, g_islands$first_locus, g_islands$last_locus))
    expect_identical(
      in_island,
      which(roh_loci$prop_roh[roh_loci$group == g] >= 0.2)
    )
  }
  expect_true(all(islands$mean_prop >= 0.2))
  expect_true(all(islands$max_prop >= islands$mean_prop))
  expect_identical(islands$start, loci$position[islands$first_locus])
  expect_identical(islands$length_bp, islands$end - islands$start + 1L)
  expect_identical(
    as.character(islands$chromosome),
    as.character(loci$chromosome[islands$last_locus])
  )
  # no islands
  none <- gt_roh_islands(sim_roh, loci, threshold = 1)
  expect_identical(nrow(none), 0L)
  expect_identical(names(none), names(islands))
  expect_error(gt_roh_islands(sim_roh, loci, threshold = 0), "threshold")
})
