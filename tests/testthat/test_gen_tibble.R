# create file
test_indiv_meta <- data.frame(
  id = c("a", "b", "c"),
  population = c("pop1", "pop1", "pop2")
)
test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, 0, 0, 0),
  c(2, 2, 0, 0, 1, 1)
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
test_that("create gen_tibble from dfs", {
  expect_true(inherits(test_gt, "gen_tbl"))
  # we can extract the genotypes correctly
  extracted_genotypes <- test_gt %>% show_genotypes()
  expect_true(all(extracted_genotypes == test_genotypes))
  # extract them from the list directly
  expect_true(all(show_genotypes(test_gt$genotypes) == test_genotypes))
  # we can extract the loci correctly
  extracted_loci <- test_gt %>% show_loci()
  # remove the index in the big file
  # expect_identical(show_loci(test_gt$genotypes) %>% select(-big_index),
  #                                                     as_tibble(test_loci))
  expect_identical(
    show_loci(test_gt$genotypes) %>%
      select(c(-big_index, -chromosome)),
    as_tibble(test_loci[, c(
      "name",
      "position",
      "genetic_dist",
      "allele_ref",
      "allele_alt"
    )])
  )
  # example of dropping the genotypes, leading to a change in class
  test_drop <- test_gt %>% select(-genotypes)
  expect_false(inherits(test_drop, "gen_tbl"))
})


# now create it directly from the dfs
test_that("create gen_tibble from dfs", {
  test_dfs_gt <- gen_tibble(
    test_genotypes,
    indiv_meta = test_indiv_meta,
    loci = test_loci,
    quiet = TRUE
  )
  # because of the different backing file info, we cannot use identical on
  # the whole object
  expect_true(identical(show_genotypes(test_gt), show_genotypes(test_dfs_gt)))
  expect_true(identical(show_loci(test_gt), show_loci(test_dfs_gt)))
  expect_true(identical(
    test_gt %>% select(-genotypes),
    test_dfs_gt %>% select(-genotypes)
  ))
})


test_that("gen_tibble catches invalid alleles", {
  test_loci_wrong <- test_loci
  test_loci_wrong$allele_alt[1] <- "N"
  expect_error(
    test_dfs_gt <- gen_tibble(
      test_genotypes,
      indiv_meta = test_indiv_meta,
      loci = test_loci_wrong,
      quiet = TRUE,
      "valid alleles are"
    )
  )
  # now add N to the valid alleles
  test_dfs_gt <- gen_tibble(
    test_genotypes,
    indiv_meta = test_indiv_meta,
    loci = test_loci_wrong,
    valid_alleles = c("A", "C", "T", "G", "N"),
    quiet = TRUE
  )
  expect_true("N" %in% show_loci(test_dfs_gt)$allele_alt)
  # but if we add to missing values it should be turned into a zero
  test_dfs_gt <- gen_tibble(
    test_genotypes,
    indiv_meta = test_indiv_meta,
    loci = test_loci_wrong,
    missing_alleles = c("0", ".", "N"),
    quiet = TRUE
  )
  expect_false("N" %in% show_loci(test_dfs_gt)$allele_alt)
  expect_true(is.na(show_loci(test_dfs_gt)$allele_alt[1]))
  # and finally throw an error if we try to use 0 as a missing value
  expect_error(
    test_dfs_gt <- gen_tibble(
      test_genotypes,
      indiv_meta = test_indiv_meta,
      loci = test_loci_wrong,
      valid_alleles = c("A", "C", "T", "G", "0"),
      quiet = TRUE
    ),
    "cannot be a valid allele"
  )
})

test_that("valid alleles check for .character method", {
  expect_error(
    pop_b <-
      gen_tibble(
        system.file("extdata/pop_b.bed", package = "tidypopgen"),
        backingfile = tempfile(),
        valid_alleles = c("A", "C", "T", "G", "0"),
        quiet = TRUE
      ),
    "cannot be a valid allele"
  )
})

test_that("error if genotypes and loci tables differ in length", {
  test_loci_short <- test_loci[1:5, ]
  expect_error(
    test_dfs_gt <- gen_tibble(
      test_genotypes,
      indiv_meta = test_indiv_meta,
      loci = test_loci_short,
      quiet = TRUE
    ),
    paste(
      "there is a mismatch between the number of loci in the genotype",
      "table x and in the loci table"
    )
  )
})

test_that("stopifnot_gen_tibble catches invalid tibbles", {
  test_gt_wrong_genotypes <- test_gt %>%
    mutate(genotypes = c(1:3))
  expect_error(
    stopifnot_gen_tibble(test_gt_wrong_genotypes),
    "the genotypes column is not of class vctrs_bigSNP"
  )

  test_gt <- test_gt %>% select(-genotypes)
  expect_error(
    indiv_missingness(test_gt),
    "'genotypes' column is missing"
  )
})

test_that("sex is stored in gen_tibble if not missing", {
  # load a gt from .bed file containing sex info
  pop_b <-
    gen_tibble(
      system.file("extdata/pop_b.bed", package = "tidypopgen"),
      backingfile = tempfile(),
      quiet = TRUE
    )
  fam <- read.table(system.file("extdata/pop_b.fam", package = "tidypopgen"),
    stringsAsFactors = FALSE
  )
  fam$V5 <- ifelse(fam$V5 == 1, "male", ifelse(fam$V5 == 2, "female", NA))
  expect_equal(as.character(pop_b$sex), as.character(fam$V5))
})

test_that("if order of loci is changed, order of genotypes also changes", {
  pop_b <-
    gen_tibble(
      system.file("extdata/pop_b.bed", package = "tidypopgen"),
      backingfile = tempfile(),
      quiet = TRUE
    )
  # original genotypes
  pop_b_gen <- show_genotypes(pop_b)

  # now scramble the loci
  set.seed(123)
  random_order <- sample(1:17)
  show_loci(pop_b) <- pop_b %>%
    select_loci(all_of(random_order)) %>%
    show_loci()

  # reorder the original genotypes according to 'random_order'
  pop_b_gen_reordered <- pop_b_gen[, random_order]

  # check that genotypes are now reordered according to random order
  expect_equal(pop_b_gen_reordered, show_genotypes(pop_b))
})

test_that("gen_tibble does not accept character matrix", {
  test_genotypes_c <- rbind(
    c("1", "1", "0", "1", "1", "0"),
    c("2", "1", "0", "0", "0", "0"),
    c("2", "2", "0", "0", "1", "1")
  )
  expect_error(
    test_dfs_gt <- gen_tibble(
      test_genotypes_c,
      indiv_meta = test_indiv_meta,
      loci = test_loci,
      quiet = TRUE
    ),
    "'x' is not a numeric matrix of integers"
  )
})

test_that("gen_tibble wrong filetype error", {
  expect_error(
    test_dfs_gt <-
      gen_tibble(system.file(
        "extdata/related/test_king.kin0",
        package = "tidypopgen"
      )),
    "x should be a valid file path pointing"
  )
})

test_that("gen_tibble loci is dataframe or tbl", {
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(3, 5, 65, 343, 23, 456)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  wrong_loci_matrix <- as.matrix(test_loci)

  expect_error(
    test_dfs_gt <- gen_tibble(
      test_genotypes,
      indiv_meta = test_indiv_meta,
      loci = wrong_loci_matrix,
      quiet = TRUE
    ),
    "loci must be a data.frame or a tibble"
  )
})

test_that("gen_tibble requires id column", {
  wrong_indiv_meta <- data.frame(
    x = c("a", "b", "c"),
    y = c("pop1", "pop1", "pop2")
  )
  expect_error(
    test_dfs_gt <- gen_tibble(
      test_genotypes,
      indiv_meta = wrong_indiv_meta,
      loci = test_loci,
      quiet = TRUE
    ),
    "indiv_meta must have the following columns: id"
  )
})

test_that("gen_tibble indiv_meta is list, dataframe, or tbl", {
  wrong_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  wrong_indiv_meta_matrix <- as.matrix(wrong_indiv_meta)

  expect_error(
    test_dfs_gt <- gen_tibble(
      test_genotypes,
      indiv_meta = wrong_indiv_meta_matrix,
      loci = test_loci,
      quiet = TRUE
    ),
    "indiv_meta must be a data.frame or a tibble"
  )
})

test_that("gen_tibble identifies wrong dimensions in genotypes", {
  wrong_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0)
  )
  expect_error(
    test_dfs_gt <- gen_tibble(
      wrong_genotypes,
      indiv_meta = test_indiv_meta,
      loci = test_loci,
      quiet = TRUE
    ),
    paste(
      "there is a mismatch between the number of individuals in the genotype",
      "table x and in the indiv_meta table"
    )
  )
})

test_that("gen_tibble identifies wrong loci table columns", {
  wrong_loci <- data.frame(
    a = paste0("rs", 1:6),
    b = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    c = as.integer(c(3, 5, 65, 343, 23, 456)),
    d = as.double(rep(0, 6)),
    e = c("A", "T", "C", "G", "C", "T"),
    f = c("T", "C", NA, "C", "G", "A")
  )
  expect_error(
    test_dfs_gt <- gen_tibble(
      test_genotypes,
      indiv_meta = test_indiv_meta,
      loci = wrong_loci,
      quiet = TRUE
    ),
    "loci must have the following columns:"
  )
})


test_that("new infrastructure for double NA", {
  # First, the case where alleles are missing
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
    allele_ref = c("A", "T", NA, "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  expect_error(gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  ), paste(
    "Some loci are missing both reference and alternate alleles.",
    "Genotypes are not missing"
  ))

  # Test with a gt from file
  expect_error(gen_tibble(
    x = test_path("testdata/plink_doubleNA", "plink_doubleNA.bed"),
    quiet = TRUE,
    backingfile = tempfile()
  ), paste(
    "Some loci are missing both reference and alternate alleles.",
    "Genotypes are not missing"
  ))


  test_genotypes <- rbind(
    c(1, 1, NA, 1, 1, 2),
    c(2, 1, NA, NA, 0, NA),
    c(2, 2, NA, 0, 1, NA)
  )
  expect_warning(test_tibble <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  ), "Your data contain loci with no genotypes or allele information")

  # Test with gt from file
  file <- tempfile()
  plink_files <- gt_as_plink(test_tibble, file = file)

  expect_warning(
    gen_tibble(
      x = plink_files, quiet = TRUE,
      backingfile = tempfile()
    ),
    "Your data contain loci with no genotypes or allele information"
  )
})

test_that("gen_tibble allow_duplicates finds duplicates", {
  # First, the case where alleles are missing
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
    name = paste0("rs", rep(1, 6)),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "G", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  expect_error(test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  ), "Your data contain duplicated locus names.")

  # But no error if allow_duplicates = TRUE
  expect_warning(test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    allow_duplicates = TRUE
  ), "Your data contain duplicated locus names.")


  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 3, 65, 343, 23, 456),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "G", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  expect_error(test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  ), "Your data contain duplicated loci.")

  # But no error if allow_duplicates = TRUE
  expect_warning(test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    allow_duplicates = TRUE
  ), "Your data contain duplicated loci.")
})


# Windows prevents the deletion of the backing file. #nolint start
# It's something to do with the memory mapping
# library used by bigsnpr
# test_that("on error, we remove the old files",{
#   # create file
#   test_indiv_meta <- data.frame (id=c("a","b","c"),
#                                  population = c("pop1","pop1","pop2"))
#   test_genotypes <- rbind(c(1,1,0,1,1,0),
#                           c(2,1,0,0,0,0),
#                           c(2,2,0,0,1,1))
#   test_loci <- data.frame(name=paste0("rs",1:6),
#                           chromosome=paste0("chr",c(1,1,1,1,2,2)),
#                           position=as.integer(c(3,5,65,343,23,456)),
#                           genetic_dist = as.double(rep(0,6)),
#                           allele_ref = c("A","T","C","G","C","T"),
#                           allele_alt = c("T","C", NA,"C","G","A"))
#   test_loci_wrong <- test_loci
#   test_loci_wrong$allele_alt[1] <- "N"
#   this_bkfile <- tempfile()
#   expect_error(test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = test_indiv_meta,
#                                          loci = test_loci_wrong,
#                                          backingfile = this_bkfile,
#                                          quiet = TRUE),"valid alleles are")
#   expect_false(file.exists(paste0(this_bkfile,".bk")))
#   test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = test_indiv_meta,
#                             loci = test_loci,
#                             backingfile = this_bkfile,
#                             quiet = TRUE)
#   expect_true(file.exists(paste0(this_bkfile,".bk")))
# }) #nolint end
