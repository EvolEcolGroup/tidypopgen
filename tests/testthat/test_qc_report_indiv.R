test_that("qc_report_indiv filters relatives within groups",{
  # Create gentibble
  bed_path <- system.file("extdata/related/families.bed", package = "tidypopgen")
  families_bigsnp_path <- bigsnpr::snp_readBed(bed_path, backingfile = tempfile())
  families <- gen_tibble(families_bigsnp_path, quiet = TRUE, valid_alleles = c("1","2"))

  ungrouped <- qc_report_indiv(families, kings_threshold = 0.2)

  # Put relatives in different populations

  #individuals 11 and 12
  #individuals 9 and 10

  families$population <- c(rep("pop1",5), rep("pop2",4),"pop1","pop2","pop1")
  # group by population to calculated relatedness within each population
  families <- families %>% group_by(population)
  # all individuals should now be retained
  grouped <- qc_report_indiv(families, kings_threshold = 0.2)
  myplots <- autoplot(grouped, type = "relatedness",kings_threshold = 0.2)
  expect_equal(length(myplots),2)
  expect_true(all(grouped$to_keep))

  families$population <- c(rep("pop1",6), rep("pop2",4),rep("pop3",2))

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
  expect_equal(as.vector(grouped$to_keep), c(pop1_report$to_keep,pop2_report$to_keep,pop3_report$to_keep))
})

test_that("qc_report_indiv$to_keep is correctly ordered",{
  # Create gentibble
  bed_path <- system.file("extdata/related/families.bed", package = "tidypopgen")
  families_bigsnp_path <- bigsnpr::snp_readBed(bed_path, backingfile = tempfile())
  families <- gen_tibble(families_bigsnp_path, quiet = TRUE, valid_alleles = c("1","2"))
  # create a population
  families$population <- c(rep("pop1",6), rep("GroupA",6))
  # make it a factor
  families$population <- as.factor(families$population)
  # GroupA is level 1, and pop1 is level 2
  expect_equal(levels(families$population), c("GroupA","pop1"))
  # But GroupA are individuals 6:12 and pop1 are individuals 1:6
  # Grouping levels are therefore not in line with order in $population
  families <- families %>% group_by(population)
  # qc_report_indiv should return the correct order
  grouped <- qc_report_indiv(families, kings_threshold = 0.2)

  families <- families %>% ungroup()
  # subset families to each population
  pop1 <- filter(families, population == "pop1")
  GroupA <- filter(families, population == "GroupA")
  # run separate qc reports
  pop1_report <- qc_report_indiv(pop1, kings_threshold = 0.2)
  GroupA_report <- qc_report_indiv(GroupA, kings_threshold = 0.2)
  # The correct individuals are removed
  expect_equal(as.vector(grouped$to_keep), c(pop1_report$to_keep,GroupA_report$to_keep))
})

