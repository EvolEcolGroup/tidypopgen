# Create gentibble
bed_path <- system.file("extdata/related/families.bed", package = "tidypopgen")
families_bigsnp_path <- bigsnpr::snp_readBed(bed_path, backingfile = tempfile())
families <- gen_tibble(
  families_bigsnp_path,
  quiet = TRUE,
  valid_alleles = c("1", "2")
)


test_that("qc_report_indiv filters relatives within groups", {
  ungrouped <- qc_report_indiv(families, kings_threshold = 0.2)

  # Put relatives in different populations

  # individuals 11 and 12
  # individuals 9 and 10

  families$population <- c(
    rep("pop1", 5),
    rep("pop2", 4),
    "pop1",
    "pop2",
    "pop1"
  )
  # group by population to calculated relatedness within each population
  families <- families %>% group_by(population)
  # all individuals should now be retained
  grouped <- qc_report_indiv(families, kings_threshold = 0.2)
  myplots <- autoplot(grouped, type = "relatedness", kings_threshold = 0.2)
  expect_equal(length(myplots), 2)
  expect_true(all(grouped$to_keep))

  families$population <- c(rep("pop1", 6), rep("pop2", 4), rep("pop3", 2))

  # subset families to each population
  pop1 <- filter(families, population == "pop1")
  pop2 <- filter(families, population == "pop2")
  pop3 <- filter(families, population == "pop3")

  # run separate qc reports
  pop1_report <- qc_report_indiv(pop1, kings_threshold = 0.2)
  pop2_report <- qc_report_indiv(pop2, kings_threshold = 0.2)
  pop3_report <- qc_report_indiv(pop3, kings_threshold = 0.2)

  families <- families %>% group_by(population)
  grouped <- qc_report_indiv(families, kings_threshold = 0.2)
  # The correct individuals are removed
  expect_equal(
    as.vector(grouped$to_keep),
    c(pop1_report$to_keep, pop2_report$to_keep, pop3_report$to_keep)
  )

  # check if there is only one individual in a population
  families[families$id == "12", "population"] <- "pop4"
  res <- qc_report_indiv(families, kings_threshold = 0.2)

  families <- families %>% ungroup()
  pop3 <- filter(families, population == "pop3")
  pop3_report <- qc_report_indiv(pop3, kings_threshold = 0.2)
  pop4 <- filter(families, population == "pop4")
  pop4_report <- qc_report_indiv(pop4, kings_threshold = 0.2)

  expect_equal(
    as.vector(res$to_keep),
    c(
      pop1_report$to_keep,
      pop2_report$to_keep,
      pop3_report$to_keep,
      pop4_report$to_keep
    )
  )
})

test_that("qc_report_indiv$to_keep is correctly ordered", {
  # create a population
  families$population <- c(rep("pop1", 6), rep("GroupA", 6))
  # make it a factor
  families$population <- as.factor(families$population)
  # GroupA is level 1, and pop1 is level 2
  expect_equal(levels(families$population), c("GroupA", "pop1"))
  # But GroupA are individuals 6:12 and pop1 are individuals 1:6
  # Grouping levels are therefore not in line with order in $population
  families <- families %>% group_by(population)
  # qc_report_indiv should return the correct order
  grouped <- qc_report_indiv(families, kings_threshold = 0.2)

  families <- families %>% ungroup()
  # subset families to each population
  pop1 <- filter(families, population == "pop1")
  group_a <- filter(families, population == "GroupA")
  # run separate qc reports
  pop1_report <- qc_report_indiv(pop1, kings_threshold = 0.2)
  group_a_report <- qc_report_indiv(group_a, kings_threshold = 0.2)
  # The correct individuals are removed
  expect_equal(
    as.vector(grouped$to_keep),
    c(pop1_report$to_keep, group_a_report$to_keep)
  )
})

test_that("autoplot qc_report_indiv list for each population", {
  ungrouped <- qc_report_indiv(families, kings_threshold = 0.2)
  expect_error(autoplot(ungrouped, type = "blah"), "'arg' should be one of")
  plot1 <- autoplot(ungrouped, type = "relatedness", kings_threshold = 0.2)
  expect_true(inherits(plot1, "gg"))
  expect_true(inherits(plot1, "ggplot"))

  families$population <- c(
    rep("pop1", 5),
    rep("pop2", 4),
    "pop1",
    "pop2",
    "pop1"
  )
  families <- families %>% group_by(population)
  grouped <- qc_report_indiv(families, kings_threshold = 0.2)

  plot2 <- autoplot(grouped, type = "relatedness", kings_threshold = 0.2)
  # expect plot2 is a list, one plot for each population
  expect_equal(names(plot2), c("pop1", "pop2"))
})

test_that("non-numeric kings_threshold arguments ", {
  expect_error(
    qc_report_indiv(families, kings_threshold = "blah"),
    "kings_threshold must be a numeric or one of"
  )

  expect_equal(
    qc_report_indiv(families, kings_threshold = "first")$to_keep,
    qc_report_indiv(families, kings_threshold = 0.177)$to_keep
  )
  expect_equal(
    qc_report_indiv(families, kings_threshold = "second")$to_keep,
    qc_report_indiv(families, kings_threshold = 0.088)$to_keep
  )


  # add population
  families$population <- c(
    rep("pop1", 5),
    rep("pop2", 4),
    "pop1",
    "pop2",
    "pop1"
  )
  # group
  families <- families %>% group_by(population)
  # test after grouping
  expect_error(
    qc_report_indiv(families, kings_threshold = "blah"),
    "kings_threshold must be a numeric or one of"
  )
  expect_equal(
    qc_report_indiv(families, kings_threshold = "first")$to_keep,
    qc_report_indiv(families, kings_threshold = 0.177)$to_keep
  )
  expect_equal(
    qc_report_indiv(families, kings_threshold = "second")$to_keep,
    qc_report_indiv(families, kings_threshold = 0.088)$to_keep
  )

  # test autoplot
  report1 <- qc_report_indiv(families, kings_threshold = "first")
  expect_equal(
    autoplot(report1,
      type = "relatedness",
      kings_threshold = "first"
    ),
    autoplot(report1,
      type = "relatedness",
      kings_threshold = 0.177
    )
  )

  report2 <- qc_report_indiv(families, kings_threshold = "second")
  expect_equal(
    autoplot(report2,
      type = "relatedness",
      kings_threshold = "second"
    ),
    autoplot(report2,
      type = "relatedness",
      kings_threshold = 0.088
    )
  )
})

test_that("qc_report_indiv for pseudohaploid data", {
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 2, 0, NA, 0, 0),
    c(2, NA, 0, 0, 1, 1),
    c(2, 0, 0, 2, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(0, 0, 0, 0, NA, 2),
    c(0, 1, 2, 0, 1, NA)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e", "f", "g"),
    population = c("pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3")
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(3, 5, 65, 343, 23, 456)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  test_gt_pseudo <- gt_pseudohaploid(test_gt)
  indiv_report <- qc_report_indiv(test_gt_pseudo)

  test_gt_pseudo <- test_gt_pseudo %>% group_by(population)
  indiv_report_grouped <- qc_report_indiv(test_gt_pseudo)

  expect_equal(indiv_report, indiv_report_grouped)

  # autoplot relatedness
  expect_error(
    autoplot(indiv_report,
      type = "relatedness",
      kings_threshold = 0.1
    ),
    "not available for pseudohaploid"
  )

  # autoplot scatter
  expect_error(
    autoplot(indiv_report,
      type = "scatter",
      kings_threshold = 0.1
    ),
    "not available for pseudohaploid"
  )

  # autoplot with kings threshold only
  expect_error(
    autoplot(indiv_report, kings_threshold = 0.177),
    "not available for pseudohaploid data"
  )

  expect_s3_class(
    autoplot(indiv_report,
      type = "histogram"
    ),
    "ggplot"
  )

  # autoplot for only pseudohaploid samples should be one pane
  test_gt_pseudo <- test_gt_pseudo %>% filter(indiv_ploidy(genotypes) == 1)
  indiv_report <- qc_report_indiv(test_gt_pseudo)
  expect_s3_class(
    autoplot(indiv_report,
      type = "histogram"
    ),
    "ggplot"
  )
})

test_that("histogram autoplot fails for non-pseudohaploid", {
  report <- qc_report_indiv(families, kings_threshold = 0.2)
  expect_error(
    autoplot(report, type = "histogram"),
    "only available for pseudohaploid data"
  )
})
