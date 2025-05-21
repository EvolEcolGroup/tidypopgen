test_that("autoplot bug with NA values", {
  # Create example
  example_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e"),
    population = c("pop1", "pop1", "pop2", "pop2", "pop2")
  )
  example_genotypes <- rbind(
    c(1, NA, 0, 1, 1, 0),
    c(2, NA, 0, 0, NA, 0),
    c(1, NA, 0, 0, 1, 1),
    c(0, NA, 0, 1, 2, 1),
    c(1, NA, NA, 2, 1, 0)
  )
  example_loci <- data.frame(
    name = c("rs1", "rs2", "rs3", "rs4", "x1", "x2"),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
    genetic_dist = c(0, 0, 0, 0, 0, 0),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  example_gt <- gen_tibble(example_genotypes,
    indiv_meta = example_indiv_meta,
    loci = example_loci,
    backingfile = tempfile(),
    quiet = TRUE
  )
  example_gt <- example_gt %>% group_by(population)

  ex_loci_report <- example_gt %>%
    qc_report_loci()

  expect_message(
    {
      pdf(NULL)
      on.exit(dev.off())
      plt <- autoplot(ex_loci_report, type = "overview")
    },
    "One or more loci are missing for every individual"
  )
  expect_false("rs2" %in% rownames(plt$New_data))
})
