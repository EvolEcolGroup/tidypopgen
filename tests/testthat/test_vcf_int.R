test_that("using integers as names works", {
  vcf_path <-
    system.file("extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
      package = "tidypopgen"
    )
  anole_gt <-
    gen_tibble(vcf_path,
      quiet = TRUE, backingfile = tempfile("anolis_"),
      names_as_int = TRUE
    )
  pops_path <- system.file("extdata/anolis/punctatus_n46_meta.csv",
    package = "tidypopgen"
  )
  pops <- read.csv(pops_path)
  anole_gt <- anole_gt %>% left_join(pops, by = "id")
  # check that names are integers
  expect_true(is.integer(show_loci(anole_gt)$name))
  anole_gt <- gt_add_sf(anole_gt, c("longitude", "latitude"))
  expect_true(is.integer(show_loci(anole_gt)$name))
  # test a pca
  anole_gt <- gt_impute_simple(anole_gt, method = "mode")
  anole_pca <- anole_gt %>% gt_pca_randomSVD(k = 30)
  expect_true(is.integer(row_names(anole_pca$v)))
  anole_pca <- anole_gt %>% gt_pca_partialSVD(k = 30)
  expect_true(is.integer(row_names(anole_pca$v)))
  # tidy the pca
  loadings_pca <- tidy(anole_pca, matrix = "loadings")
  # check that we have a an integer column for "column"
  expect_true(is.integer(loadings_pca$column))
})

# helper to create a gen_tibble with integer loci names that are not 1:n
# (to make sure that names are not confused with positions)
make_int_gt <- function(names_as_int = TRUE) {
  vcf_path <-
    system.file("extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
      package = "tidypopgen"
    )
  anole_gt <- gen_tibble(vcf_path,
    quiet = TRUE, backingfile = tempfile("anolis_"),
    names_as_int = names_as_int
  )
  anole_gt <- anole_gt %>% select_loci(seq(2, count_loci(anole_gt), by = 2))
  anole_gt <- gt_update_backingfile(anole_gt,
    backingfile = tempfile("anolis_sub_"), quiet = TRUE
  )
  anole_gt$population <- rep(c("a", "b", "c"), length.out = nrow(anole_gt))
  anole_gt
}

for (names_as_int in c(TRUE, FALSE)) {
  test_that(paste(
    "grouped loci matrices keep loci names, names_as_int =",
    names_as_int
  ), {
    anole_gt <- make_int_gt(names_as_int)
    expect_equal(is.integer(loci_names(anole_gt)), names_as_int)
    anole_grp <- anole_gt %>% group_by(population)
    loci_funs <- list(loci_alt_freq, loci_maf, loci_missingness, loci_hwe)
    for (this_fun in loci_funs) {
      res <- this_fun(anole_grp, type = "matrix")
      expect_identical(row_names(res), loci_names(anole_gt))
      expect_identical(col_names(res), c("a", "b", "c"))
    }
    res <- loci_pi(anole_grp, type = "matrix")
    expect_identical(row_names(res), loci_names(anole_gt))
    # maf values are consistent with the ungrouped estimates
    maf_mat <- loci_maf(anole_grp, type = "matrix")
    expect_equal(
      as.vector(maf_mat[, "b"]),
      loci_maf(anole_gt %>% filter(population == "b"))
    )
  })

  test_that(paste(
    "nwise_pop_pbs keeps loci names, names_as_int =",
    names_as_int
  ), {
    anole_gt <- make_int_gt(names_as_int)
    anole_grp <- anole_gt %>% group_by(population)
    pbs_mat <- nwise_pop_pbs(anole_grp, type = "matrix", return_fst = TRUE)
    expect_identical(row_names(pbs_mat), loci_names(anole_gt))
    expect_true("fst_a.b" %in% colnames(pbs_mat))
    pbs_tidy <- nwise_pop_pbs(anole_grp, type = "tidy", return_fst = TRUE)
    expect_identical(unique(pbs_tidy$loci), loci_names(anole_gt))
    expect_equal(
      pbs_tidy$value[pbs_tidy$stat_name == "pbs_a.b.c"],
      as.vector(pbs_mat[, "pbs_a.b.c"])
    )
    expect_equal(
      pbs_tidy$value[pbs_tidy$stat_name == "fst_a.b"],
      as.vector(pbs_mat[, "fst_a.b"])
    )
  })

  test_that(paste(
    "pca and dapc keep loci names, names_as_int =",
    names_as_int
  ), {
    anole_gt <- make_int_gt(names_as_int)
    anole_gt <- gt_impute_simple(anole_gt, method = "mode")
    anole_pca <- gt_pca_partialSVD(anole_gt, k = 5)
    expect_identical(row_names(anole_pca$v), loci_names(anole_gt))
    # tidy loadings
    loadings_pca <- tidy(anole_pca, matrix = "loadings")
    expect_identical(unique(loadings_pca$column), loci_names(anole_gt))
    expect_equal(
      loadings_pca$value[loadings_pca$PC == 2],
      as.vector(anole_pca$v[, 2])
    )
    # augment_loci
    aug_loci <- augment_loci(anole_pca, data = anole_gt)
    expect_identical(aug_loci$name, loci_names(anole_gt))
    expect_equal(aug_loci$.loadingPC1, as.vector(anole_pca$v[, 1]))
    aug_loci <- augment_loci(anole_pca)
    expect_identical(aug_loci$.rownames, loci_names(anole_gt))
    # autoSVD only keeps a subset of loci
    anole_auto <- gt_pca_autoSVD(anole_gt,
      k = 5, roll_size = 0, verbose = FALSE
    )
    expect_identical(row_names(anole_auto$v), anole_auto$loci$name)
    # dapc
    anole_dapc <- gt_dapc(anole_pca,
      pop = as.factor(anole_gt$population),
      n_pca = 5, n_da = 2
    )
    expect_identical(row_names(anole_dapc$var.load), loci_names(anole_gt))
    expect_identical(row_names(anole_dapc$var.contr), loci_names(anole_gt))
    loadings_dapc <- tidy(anole_dapc, matrix = "loadings")
    expect_identical(unique(loadings_dapc$column), loci_names(anole_gt))
    expect_equal(
      loadings_dapc$value[loadings_dapc$LD == "LD1"],
      as.vector(anole_dapc$var.load[, 1])
    )
  })
}
