# Compare Tajima's D with missing data to pixy (Bailey et al. 2025)
# The vcf and the pixy results were generated with pixy version 2.2.4
# See data-raw/create_pixy_test_data.R for the script

# read the test data into a grouped gen_tibble
pixy_gt <- gen_tibble(
  test_path("testdata/pixy", "pixy_test.vcf"),
  backingfile = tempfile(),
  quiet = TRUE
)
pixy_pops <- utils::read.table(
  test_path("testdata/pixy", "pixy_populations.txt"),
  col.names = c("id", "population")
)
pixy_gt$population <- pixy_pops$population[match(pixy_gt$id, pixy_pops$id)]
pixy_gt <- pixy_gt %>% dplyr::group_by(population)

test_that("windows_pop_tajimas_d matches pixy", {
  pixy_windows <- utils::read.delim(
    test_path("testdata/pixy", "pixy_10000_tajima_d.txt")
  )
  # windows of 10kb with no overlap match the windows used by pixy
  window_taj <- windows_pop_tajimas_d(
    pixy_gt,
    type = "tidy",
    window_size = 10000,
    step_size = 10000,
    size_unit = "bp"
  )
  window_taj <- merge(pixy_windows, window_taj,
    by.x = c("pop", "chromosome", "window_pos_1", "window_pos_2"),
    by.y = c("group", "chromosome", "start", "end")
  )
  # all windows estimated by pixy are found
  expect_equal(nrow(window_taj), nrow(pixy_windows))
  expect_equal(window_taj$stat, window_taj$tajima_d)
  # pixy counts the number of sites with at least one genotype
  expect_equal(window_taj$n_loci, window_taj$no_sites)
})

test_that("pop_tajimas_d and jackknife_pop_tajimas_d match pixy", {
  # pixy computes the statistics for each chromosome, so we use a window
  # larger than the chromosomes
  pixy_chrom <- utils::read.delim(
    test_path("testdata/pixy", "pixy_100000_tajima_d.txt")
  )
  for (chrom in c("chr1", "chr2")) {
    pixy_gt_chrom <- pixy_gt %>%
      select_loci_if(loci_chromosomes(genotypes) == chrom)
    pixy_d <- pixy_chrom$tajima_d[pixy_chrom$chromosome == chrom]
    names(pixy_d) <- pixy_chrom$pop[pixy_chrom$chromosome == chrom]
    # genome-wide estimate
    taj_chrom <- pop_tajimas_d(pixy_gt_chrom)
    expect_equal(unlist(taj_chrom), unname(pixy_d[c("popA", "popB", "popC")]))
    # the jackknife estimate on all loci is the same
    jack_chrom <- jackknife_pop_tajimas_d(pixy_gt_chrom,
      type = "tidy",
      block_size = 10000,
      size_unit = "bp"
    )
    expect_equal(jack_chrom$estimate, unname(pixy_d[jack_chrom$group]))
  }
})
