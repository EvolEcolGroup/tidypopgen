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

})
