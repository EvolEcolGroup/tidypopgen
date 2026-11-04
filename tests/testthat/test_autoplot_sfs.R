# Tests for the autoplot() methods of site frequency spectra

sim <- ref_sfs_simulate(c(popA = 6, popB = 5, popC = 4),
  n_loci = 200,
  seed = 51
)

# the data of the first layer of a plot, after building it
layer_data_of <- function(p, i = 1) {
  ggplot2::layer_data(p, i)
}

test_that("pop_sfs tidy outputs are plotted as bar plots", {
  sfs <- pop_sfs(sim$gt, n_proj = c(8, 6, 4))
  expect_s3_class(sfs, "pop_sfs_tbl")
  p <- autoplot(sfs)
  expect_true(inherits(p, "ggplot"))
  expect_s3_class(p$facet, "FacetWrap")
  expect_identical(p$labels$x, "Minor allele count")
  # masked cells (corners and folded out cells) are not shown
  expect_identical(nrow(p$data), sum(!sfs$masked))
  expect_equal(p$data$n_sites, sfs$n_sites[!sfs$masked])
  expect_identical(
    levels(p$data$panel),
    c("popA", "popB", "popC")
  )
  # unless asked for
  p_masked <- autoplot(sfs, masked = TRUE)
  expect_identical(nrow(p_masked$data), nrow(sfs))
  # unfolded spectra count derived alleles
  p_unfolded <- autoplot(pop_sfs(sim$gt, n_proj = 4, folded = FALSE))
  expect_identical(p_unfolded$labels$x, "Derived allele count")
  # a log scale drops empty cells
  p_log <- autoplot(sfs, log = TRUE, masked = TRUE)
  expect_true(all(p_log$data$n_sites > 0))
  expect_no_error(ggplot2::ggplot_build(p_log))
  # ungrouped data give a single panel
  p_all <- autoplot(pop_sfs(dplyr::ungroup(sim$gt), n_proj = 10))
  expect_s3_class(p_all$facet, "FacetNull")
})

test_that("pairwise_pop_sfs tidy outputs are plotted as heatmaps", {
  sfs <- pairwise_pop_sfs(sim$gt, n_proj = c(8, 6, 4), folded = FALSE)
  expect_s3_class(sfs, "pairwise_sfs_tbl")
  p <- autoplot(sfs)
  expect_true(inherits(p, "ggplot"))
  expect_identical(
    levels(p$data$panel),
    c("popA (x) - popB (y)", "popA (x) - popC (y)", "popB (x) - popC (y)")
  )
  # masked and empty cells are left blank
  expect_true(all(is.na(p$data$n_sites[sfs$masked])))
  expect_true(all(is.na(p$data$n_sites[sfs$n_sites == 0])))
  # the colour scale starts at a thousandth of the largest number of sites
  shown <- sfs$n_sites[!sfs$masked & sfs$n_sites > 0]
  min_sites <- max(shown) / 1000
  expect_equal(min(p$data$n_sites, na.rm = TRUE), min_sites)
  expect_equal(
    p$data$n_sites[!is.na(p$data$n_sites)],
    pmax(shown, min_sites)
  )
  # or at a given number of sites
  p_min <- autoplot(sfs, min_sites = 0.5)
  expect_equal(min(p_min$data$n_sites, na.rm = TRUE), 0.5)
  expect_error(autoplot(sfs, min_sites = -1), "positive number")
  # on a linear scale, empty cells are shown
  p_linear <- autoplot(sfs, log = FALSE)
  expect_true(any(p_linear$data$n_sites == 0, na.rm = TRUE))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("sfs objects are plotted", {
  sfs_1d <- pop_sfs(sim$gt, n_proj = c(8, 6, 4), type = "array")$popB
  p_1d <- autoplot(sfs_1d)
  expect_true(inherits(p_1d, "ggplot"))
  expect_identical(p_1d$labels$title, "popB")
  expect_equal(p_1d$data$n_sites, as.vector(sfs_1d)[2:4])
  sfs_2d <- pairwise_pop_sfs(sim$gt,
    n_proj = c(8, 6, 4), type = "array"
  )$popA.popC
  p_2d <- autoplot(sfs_2d)
  expect_true(inherits(p_2d, "ggplot"))
  expect_identical(p_2d$labels$x, "Minor allele count (popA)")
  expect_identical(p_2d$labels$y, "Minor allele count (popC)")
  expect_no_error(ggplot2::ggplot_build(p_2d))
  sfs_3d <- nwise_pop_sfs(sim$gt, n_proj = c(4, 4, 2), type = "array")
  expect_error(autoplot(sfs_3d), "one or two populations")
})

test_that("pop_sfs_preview outputs are plotted", {
  preview <- pop_sfs_preview(sim$gt)
  expect_s3_class(preview, "sfs_preview_tbl")
  p <- autoplot(preview)
  expect_true(inherits(p, "ggplot"))
  # the line of segregating sites, the best projection, and its label
  expect_identical(nrow(layer_data_of(p, 1)), nrow(preview))
  expect_identical(nrow(layer_data_of(p, 2)), 3L)
  expect_identical(
    layer_data_of(p, 3)$label,
    paste0("n_proj = ", preview$n_proj[preview$best])
  )
  # the class is kept by dplyr verbs
  preview_ab <- preview %>% dplyr::filter(group != "popC")
  expect_s3_class(preview_ab, "sfs_preview_tbl")
  expect_true(inherits(autoplot(preview_ab), "ggplot"))
})

test_that("autoplot checks its input", {
  sfs <- pop_sfs(sim$gt, n_proj = 4)
  expect_error(
    autoplot(dplyr::select(sfs, -"masked")),
    "should be the output of pop_sfs()"
  )
  expect_error(autoplot(sfs, masked = NA), "`masked` should be TRUE or FALSE")
  expect_error(autoplot(sfs, colour = "red"), "must be empty")
})
