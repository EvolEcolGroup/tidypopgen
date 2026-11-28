# the fixture (data-raw/validate/create_plink_sexcheck_test_data.R): 18
# individuals with loci on 1, X, Y and MT, and the outputs of PLINK 1.9
fixture <- test_path("testdata", "plink_sexcheck", "sexcheck")

# the sex codes of PLINK (1 male, 2 female, 0 unknown) as a factor
plink_sex <- function(codes) {
  factor(c("male", "female")[ifelse(codes == 0, NA, codes)],
    levels = c("female", "male")
  )
}

test_that("indiv_sex_check matches PLINK 1.9 --check-sex", {
  sex_gt <- gen_tibble(paste0(fixture, ".bed"),
    backingfile = tempfile(), quiet = TRUE
  )
  # the daughter's parents are known (so she is not a founder)
  expect_identical(sex_gt$paternal_ID[17], "ind9")
  plink_res <- utils::read.table(paste0(fixture, ".sexcheck"), header = TRUE)
  res <- indiv_sex_check(sex_gt)
  expect_identical(res$id, plink_res$IID)
  expect_identical(res$sex, plink_sex(plink_res$PEDSEX))
  expect_identical(res$snp_sex, plink_sex(plink_res$SNPSEX))
  expect_identical(res$status, plink_res$STATUS)
  # (PLINK writes F with 4 significant digits)
  expect_equal(res$f, plink_res$F, tolerance = 1e-3)
  expect_false("y_count" %in% names(res))
  # with the counts of genotypes on Y (males need at least 5)
  plink_y <- utils::read.table(paste0(fixture, "_ycount.sexcheck"),
    header = TRUE
  )
  res_y <- indiv_sex_check(sex_gt, ycount = TRUE, min_male_ycount = 5)
  expect_identical(res_y$snp_sex, plink_sex(plink_y$SNPSEX))
  expect_identical(res_y$status, plink_y$STATUS)
  expect_identical(res_y$y_count, plink_y$YCOUNT)
  # groups are ignored
  expect_identical(
    indiv_sex_check(dplyr::group_by(sex_gt, sex)),
    res
  )
})

test_that("indiv_sex_check needs a polymorphic X, and checks its inputs", {
  sex_gt <- gen_tibble(paste0(fixture, ".bed"),
    backingfile = tempfile(), quiet = TRUE
  )
  # without X
  expect_error(
    sex_gt %>%
      select_loci_if(loci_chromosome_types(genotypes) != "X") %>%
      indiv_sex_check(),
    "no loci on the X"
  )
  # with a monomorphic X (the first locus on X)
  x_loci <- which(loci_chromosome_types(sex_gt) == "X")
  expect_error(
    sex_gt %>%
      select_loci(c(1, x_loci[1])) %>%
      indiv_sex_check(),
    "at least one polymorphic locus"
  )
  # recorded sexes as characters or PLINK's codes
  sex_gt$sex <- c("M", "F")[as.integer(sex_gt$sex == "female") + 1]
  expect_identical(
    indiv_sex_check(sex_gt)$sex[1:3],
    plink_sex(c(2, 2, 1))
  )
  sex_gt$sex <- c(2, 2, 1, rep(0, 15))
  expect_identical(
    indiv_sex_check(sex_gt)$sex[1:4],
    plink_sex(c(2, 2, 1, 0))
  )
  sex_gt$sex <- "hermaphrodite"
  expect_error(indiv_sex_check(sex_gt), "should give 'female' or 'male'")
})

test_that("indiv_sex_check uses Z in ZW systems (females are hemizygous)", {
  sex_gt <- gen_tibble(paste0(fixture, ".bed"),
    backingfile = tempfile(), quiet = TRUE
  )
  # the X becomes a Z (and Y a W), and the sexes are swapped
  swap_sex <- function(x) {
    factor(unname(c(female = "male", male = "female")[as.character(x)]),
      levels = c("female", "male")
    )
  }
  zw_gt <- gt_set_chromosome_types(sex_gt, types = c(X = "Z", Y = "W"))
  zw_gt$sex <- swap_sex(sex_gt$sex)
  res <- indiv_sex_check(sex_gt)
  res_zw <- indiv_sex_check(zw_gt)
  expect_identical(res_zw$f, res$f)
  expect_identical(res_zw$status, res$status)
  expect_identical(res_zw$snp_sex, swap_sex(res$snp_sex))
})

test_that("haploid genotypes of a VCF are read as PLINK reads them", {
  # males are haploid on X, Y and MT in the VCF, which gives the same
  # genotypes as PLINK's .bed of that VCF (haploid genotypes as homozygotes)
  vcf_gt <- gen_tibble(paste0(fixture, ".vcf"),
    backingfile = tempfile(), quiet = TRUE
  )
  bed_gt <- gen_tibble(paste0(fixture, "_vcf.bed"),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(indiv_ploidy(vcf_gt), rep(2L, 18))
  expect_identical(show_genotypes(vcf_gt), show_genotypes(bed_gt))
})
