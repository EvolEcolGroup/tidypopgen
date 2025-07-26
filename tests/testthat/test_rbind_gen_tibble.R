raw_path_pop_b <- system.file("extdata/pop_b.bed", package = "tidypopgen")
bigsnp_path_b <- bigsnpr::snp_readBed(
  raw_path_pop_b,
  backingfile = tempfile("test_b_")
)
pop_b_gt <- gen_tibble(bigsnp_path_b, quiet = TRUE)
# target file
raw_path_pop_a <- system.file("extdata/pop_a.bed", package = "tidypopgen")
bigsnp_path_a <- bigsnpr::snp_readBed(
  raw_path_pop_a,
  backingfile = tempfile("test_a_")
)
pop_a_gt <- gen_tibble(bigsnp_path_a, quiet = TRUE)
# #create merge
merged_gen <- rbind.gen_tbl(
  pop_b_gt,
  pop_a_gt,
  flip_strand = TRUE,
  quiet = TRUE,
  backingfile = tempfile()
)

test_that("merge combines datasets correctly", {
  genotypes <- show_genotypes(merged_gen)

  # Genotypes before merging
  pop_b_geno <- show_genotypes(pop_b_gt)
  pop_a_geno <- show_genotypes(pop_a_gt)

  # Check pop_b
  pop_b_merged <- merged_gen %>%
    filter(population == "pop_b") %>%
    show_genotypes()

  # All pop_b genotypes should stay the same, as they are the reference
  expect_equal(
    pop_b_merged,
    pop_b_geno[, c(1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17)]
  )

  # Check pop_b
  pop_a_merged <- merged_gen %>%
    filter(population == "pop_a") %>%
    show_genotypes()

  # Check genotypes of non-swapped alleles are the same before and after merge
  expect_equal(pop_a_merged[, c(1, 2, 3)], pop_a_geno[, c(2, 3, 4)])

  # Check swapped alleles have genotypes swapped correctly
  expect_true(all(pop_a_merged[, 4] == c(2, 2, 0, 2, 1))) # rs11240777
  expect_true(all(pop_a_merged[, 10] == c(1, 2, 1, 1, 1))) # rs10106770
  expect_true(all(pop_a_merged[, 11] == c(2, 1, 2, 1, 2))) # rs11942835

  # Check ambiguous SNPs are dropped (remove_ambiguous TRUE by default)
  expect_false("rs1240719" %in% loci_names(merged_gen))
  expect_false("rs307354" %in% loci_names(merged_gen))
})

test_that("warning when no snps overlap", {
  # Create two datasets

  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 2),
    c(2, 1, 0, NA, 0, NA),
    c(2, 2, 0, 0, 1, NA)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
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

  test_indiv_meta2 <- data.frame(
    id = c("A", "B", "C"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes2 <- rbind(
    c(1, 1, 0, 1, 1, 2),
    c(2, 1, 0, NA, 0, NA),
    c(2, 2, 0, 0, 1, NA)
  )
  test_loci2 <- data.frame(
    name = paste0("rs", 7:12),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt2 <- gen_tibble(
    x = test_genotypes2,
    loci = test_loci2,
    indiv_meta = test_indiv_meta2,
    quiet = TRUE
  )

  report1 <- rbind_dry_run(test_gt, test_gt2, flip_strand = TRUE, quiet = TRUE)

  # merge
  expect_error(
    rbind(
      test_gt,
      test_gt2,
      flip_strand = TRUE,
      quiet = TRUE,
      backingfile = tempfile()
    ),
    "there are no loci in common between"
  )

  # Now try with the same snp ID but different CHR

  # These should merge, with CHR and POS in merged GT being equal to test_gt2

  test_indiv_meta2 <- data.frame(
    id = c("A", "B", "C"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes2 <- rbind(
    c(1, 1, 0, 1, 1, 2),
    c(2, 1, 0, NA, 0, NA),
    c(2, 2, 0, 0, 1, NA)
  )
  test_loci2 <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = c(3, 3, 3, 3, 4, 4),
    position = c(3, 5, 65, 343, 23, 456),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt2 <- gen_tibble(
    x = test_genotypes2,
    loci = test_loci2,
    indiv_meta = test_indiv_meta2,
    quiet = TRUE
  )

  report2 <- rbind_dry_run(test_gt, test_gt2, flip_strand = TRUE, quiet = TRUE)

  # merge
  merged_gt <- rbind(
    test_gt,
    test_gt2,
    flip_strand = TRUE,
    backingfile = tempfile(),
    quiet = TRUE
  )

  expect_true(all(show_loci(merged_gt)$"chromosome" == c(1, 1)))
  expect_true(all(show_loci(merged_gt)$"position" == c(5, 65)))
})


test_that("merge by position works correctly", {
  pop_a_renamed_gt <- pop_a_gt
  show_loci(pop_a_renamed_gt)$name <- paste(
    "new_name",
    1:count_loci(pop_a_renamed_gt),
    sep = "_"
  )
  # if we bind by name we should get zero (or an error)
  expect_error(
    rbind(pop_b_gt, pop_a_renamed_gt, quiet = TRUE),
    "there are no loci in common"
  )
  pos_merge <- rbind(
    pop_b_gt,
    pop_a_renamed_gt,
    flip_strand = TRUE,
    use_position = TRUE,
    quiet = TRUE
  )
  expect_true(all.equal(show_genotypes(merged_gen), show_genotypes(pos_merge)))
})

test_that("original gen_tibble objects remain the same after merge", {
  # reference

  raw_path_pop_b <- system.file("extdata/pop_b.bed", package = "tidypopgen")
  bigsnp_path_b <- bigsnpr::snp_readBed(
    raw_path_pop_b,
    backingfile = tempfile("test_b_")
  )
  pop_b_gt <- gen_tibble(bigsnp_path_b, quiet = TRUE)

  # hash before merge

  md5_before_b_rds <- tools::md5sum(bigsnp_path_b)
  bk_path_b <- sub("\\.rds$", ".bk", bigsnp_path_b)
  md5_before_b_bk <- tools::md5sum(bk_path_b)

  # copy object

  pop_b_gt_copy <- pop_b_gt

  # target

  raw_path_pop_a <- system.file("extdata/pop_a.bed", package = "tidypopgen")
  bigsnp_path_a <- bigsnpr::snp_readBed(
    raw_path_pop_a,
    backingfile = tempfile("test_a_")
  )
  pop_a_gt <- gen_tibble(bigsnp_path_a, quiet = TRUE)

  # hash before merge

  md5_before_a_rds <- tools::md5sum(bigsnp_path_a)
  bk_path_a <- sub("\\.rds$", ".bk", bigsnp_path_a)
  md5_before_a_bk <- tools::md5sum(bk_path_a)

  # copy object

  pop_a_gt_copy <- pop_a_gt

  # merge

  merged_gen <- rbind.gen_tbl(
    pop_b_gt,
    pop_a_gt,
    flip_strand = TRUE,
    quiet = TRUE,
    backingfile = tempfile()
  )

  # check objects

  expect_equal(pop_a_gt_copy, pop_a_gt)
  expect_equal(pop_b_gt_copy, pop_b_gt)

  # check attributes
  expect_equal(attributes(pop_a_gt_copy), attributes(pop_a_gt))
  expect_equal(attributes(pop_b_gt_copy), attributes(pop_b_gt))

  # hash after merge

  md5_after_b_rds <- tools::md5sum(bigsnp_path_b)
  md5_after_a_rds <- tools::md5sum(bigsnp_path_a)

  # hash after merge

  md5_after_b_bk <- tools::md5sum(bk_path_b)
  md5_after_a_bk <- tools::md5sum(bk_path_a)

  expect_equal(md5_before_a_rds, md5_after_a_rds)
  expect_equal(md5_before_b_rds, md5_after_b_rds)

  expect_equal(md5_before_a_bk, md5_after_a_bk)
  expect_equal(md5_before_b_bk, md5_after_b_bk)
})

test_that("rbind fails if gen_tibbles contain same ID", {
  pop_a_gt$id[1] <- "A"
  pop_b_gt$id[1] <- "A"

  expect_error(merged_gen <- rbind.gen_tbl(
    pop_b_gt,
    pop_a_gt,
    flip_strand = TRUE,
    quiet = TRUE,
    backingfile = tempfile()
  ), "at least one individual with the same ID")
})

test_that("rbind warning when id is duplicated in bigsnp object", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c")
  )
  test_genotypes <- rbind(
    c(2, 2, 2, 2, 2, 2),
    c(2, 2, 2, 2, 2, 2),
    c(2, 2, 2, 2, 2, 2)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt1 <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  # change genotypes
  test_genotypes2 <- rbind(
    c(0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0)
  )
  # create a second test gt with identical id
  test_gt2 <- gen_tibble(
    x = test_genotypes2,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  # add population to both tests_gt1 and test_gt2
  test_gt1 <- test_gt1 %>% dplyr::mutate(population = "pop1")
  test_gt2 <- test_gt2 %>% dplyr::mutate(population = "pop2")
  # update id according to population
  test_gt1$id <- paste(test_gt1$id, test_gt1$population, sep = "_")
  test_gt2$id <- paste(test_gt2$id, test_gt2$population, sep = "_")
  # merge
  expect_error(
    rbind(test_gt1, test_gt2),
    "The two bigsnp objects contain at least one individual "
  )
})
