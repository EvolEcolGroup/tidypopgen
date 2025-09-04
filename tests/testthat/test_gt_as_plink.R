# create file
test_indiv_meta <- data.frame(
  id = c("a", "b", "c"),
  population = c("pop1", "pop1", "pop2")
)
test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, 0, 0, 0),
  c(2, 2, 0, 0, NA, 1)
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


# this also tests show_genotypes and show_loci
test_that("write a bed file", {
  bed_path <- gt_as_plink(test_gt, file = paste0(tempfile(), ".bed"))
  # now read the file back in
  test_gt2 <- gen_tibble(bed_path, quiet = TRUE)
  ## continue here

  # because of the different backing file info
  # we cannot use identical on the whole object
  expect_true(identical(show_genotypes(test_gt), show_genotypes(test_gt2)))
  expect_true(identical(
    test_gt %>% select(-genotypes),
    test_gt2 %>% select(-genotypes)
  ))

  # check gt_as_plink converts the NA missing allele to 0
  expect_true(is.na(show_loci(test_gt2)$allele_alt[3]))

  # now write it as a ped
  ped_path <- gt_as_plink(
    test_gt,
    file = paste0(tempfile(), ".ped"),
    type = "ped"
  )
  test_gt3 <- gen_tibble(ped_path, quiet = TRUE)
  # the gen tibble from the bed and ped should contain the same information
  expect_true(all.equal(
    show_loci(test_gt3),
    show_loci(test_gt2),
    check.attributes = FALSE
  ))
  expect_true(all.equal(show_genotypes(test_gt3), show_genotypes(test_gt2)))

  # write it as raw
  raw_path <- gt_as_plink(test_gt, file = tempfile(), type = "raw")
  raw_file_test <- read.table(raw_path, header = TRUE)
  mat <- as.matrix(raw_file_test[, 7:ncol(raw_file_test)])
  mat <- unname(mat)
  expect_true(all.equal(mat, show_genotypes(test_gt)))
})

test_that("test for overwriting files", {
  temp <- tempfile()

  # write a bed
  bed_path <- gt_as_plink(test_gt, file = paste0(temp, ".bed"))

  # files exist
  expect_true(file.exists(paste0(temp, ".bed")))
  expect_true(file.exists(paste0(temp, ".bim")))
  expect_true(file.exists(paste0(temp, ".fam")))

  # now write it as a ped
  ped_path <- gt_as_plink(test_gt, file = paste0(temp, ".ped"), type = "ped")

  # files exist
  expect_true(file.exists(paste0(temp, ".bed")))
  expect_true(file.exists(paste0(temp, ".bim")))
  expect_true(file.exists(paste0(temp, ".fam")))
  expect_true(file.exists(paste0(temp, ".ped")))
  expect_true(file.exists(paste0(temp, ".map")))

  # if .ped is created first

  temp2 <- tempfile()

  # write a ped
  ped_path <- gt_as_plink(test_gt, file = paste0(temp2, ".ped"), type = "ped")
  # files exist
  expect_true(file.exists(paste0(temp, ".ped")))
  expect_true(file.exists(paste0(temp, ".map")))

  # then write is as a bed
  bed_path <- gt_as_plink(test_gt, file = paste0(temp2, ".bed"))
  # files exist
  expect_true(file.exists(paste0(temp2, ".ped")))
  expect_true(file.exists(paste0(temp2, ".map")))
  expect_true(file.exists(paste0(temp2, ".bed")))
  expect_true(file.exists(paste0(temp2, ".bim")))
  expect_true(file.exists(paste0(temp2, ".fam")))
})

test_that("gt_as_plink uses loci and indiv information from the gen_tibble", {
  test_gt$population <- c("population1", "population1", "population2")
  test_gt$id <- c("indiv1", "indiv2", "indiv3")

  loci <- show_loci(test_gt)
  loci$name <- paste0("rs_id", 1:6)
  show_loci(test_gt) <- loci

  # save using gt_as_plink
  bed_path <- gt_as_plink(test_gt, file = paste0(tempfile(), ".bed"))

  # read the bed file
  test_gt2 <- gen_tibble(bed_path, quiet = TRUE)

  expect_equal(show_loci(test_gt), show_loci(test_gt2))
  expect_equal(test_gt$id, test_gt2$id)

  # now check the same with a grouped gen_tibble
  test_gt_grouped <- test_gt %>% group_by(population)

  # write the bed
  bed_path2 <- gt_as_plink(test_gt_grouped, file = paste0(tempfile(), ".bed"))

  # read the bed
  test_gt3 <- gen_tibble(bed_path2, quiet = TRUE)

  expect_equal(show_loci(test_gt), show_loci(test_gt3))
  expect_equal(test_gt$id, test_gt3$id)
  expect_equal(test_gt$population, test_gt3$population)
})

test_that("gt_as_plink can use chr_int", {
  # save using gt_as_plink with chromosomes_as_int TRUE
  bed_path <- gt_as_plink(
    test_gt,
    file = paste0(tempfile(), ".bed"),
    chromosomes_as_int = TRUE
  )

  # read the bed file
  test_gt_chr_int <- gen_tibble(bed_path, quiet = TRUE)

  # take original bed, remove the chr
  result <- sub("^chr", "", show_loci(test_gt)$chromosome)

  # compare
  expect_equal(result, show_loci(test_gt_chr_int)$chromosome)
})

if (rlang::is_installed("vcfR")) {
  test_that("family.ID equals sample.ID from vcf", {
    # If the gen_tibble has been read in from vcf format,
    # family.ID in the resulting plink files will be the same as sample.ID.

    ####  With vcfr parser
    vcf_path <- system.file("extdata/pop_b.vcf", package = "tidypopgen")
    pop_b_vcf_gt <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      backingfile = tempfile(),
      parser = "vcfR"
    )
    # write vcf_path using gt_as_plink
    pop_b_bed <- gt_as_plink(pop_b_vcf_gt, tempfile())
    # substitute ".bed" for ".fam"
    fam_path <- gsub(".bed", ".fam", pop_b_bed)
    # read in the .fam file
    fam <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
    # check that family.ID is the same as sample.ID
    expect_true(all(fam$V1 == pop_b_vcf_gt$id))
    expect_true(all(fam$V1 == fam$V2))

    ####  With cpp parser
    vcf_path <- system.file("extdata/pop_b.vcf", package = "tidypopgen")
    pop_b_vcf_gt <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      backingfile = tempfile(),
      parser = "cpp"
    )
    # write vcf_path using gt_as_plink
    pop_b_bed <- gt_as_plink(pop_b_vcf_gt, tempfile())
    # substitute ".bed" for ".fam"
    fam_path <- gsub(".bed", ".fam", pop_b_bed)
    # read in the .fam file
    fam <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
    # check that family.ID is the same as sample.ID
    expect_true(all(fam$V1 == pop_b_vcf_gt$id))
    expect_true(all(fam$V1 == fam$V2))
  })
}

test_that("handling of duplicated loci", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, NA, 1)
  )
  # create a gt with duplicated positions
  test_loci <- data.frame(
    name = paste0("rs", c(1, 2, 1, 3:5)),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(65, 5, 65, 343, 23, 456)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("C", "T", "C", "G", "C", "T"),
    allele_alt = c(NA, "C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  # duplicated positions are registered by is_loci_table_ordered
  expect_false(is_loci_table_ordered(test_gt, ignore_genetic_dist = FALSE))
  expect_error(
    is_loci_table_ordered(
      test_gt,
      ignore_genetic_dist = FALSE,
      error_on_false = TRUE
    ),
    "Your loci are not sorted within chromosomes"
  )
  # but currently we can still write the plink fileset and read back
  # into a gt without a warning
  file <- gt_as_plink(
    test_gt,
    file = paste0(tempfile(), ".bed"),
    chromosomes_as_int = TRUE
  )
  test_gt_reload <- gen_tibble(file, backingfile = tempfile(), quiet = TRUE)

  # create a gt with adjacent duplicated positions
  test_loci <- data.frame(
    name = paste0("rs", c(1, 2, 2, 3:5)),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(3, 65, 65, 343, 23, 456)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "C", "C", "G", "C", "T"),
    allele_alt = c("T", NA, NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  # duplicated positions are registered by is_loci_table_ordered
  expect_false(is_loci_table_ordered(test_gt, ignore_genetic_dist = FALSE))
  # but currently we can still write the plink fileset and read back
  # into a gt without a warning
  file <- gt_as_plink(
    test_gt,
    file = paste0(tempfile(), ".bed"),
    chromosomes_as_int = TRUE
  )
  test_gt_reload <- gen_tibble(file, backingfile = tempfile(), quiet = TRUE)
})
