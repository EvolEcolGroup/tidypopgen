test_that("jackknife_pop_tajimas_d leave-one-out matches subsets", {
  test_gt <- make_jack_test_gt()
  block_id <- .jackknife_block_ids(show_loci(test_gt),
    block_size = 1, size_unit = "cM"
  )
  jack_list <- jackknife_pop_tajimas_d(test_gt,
    type = "list", block_size = 1, size_unit = "cM"
  )
  # estimate is the genome-wide tajima's d
  expect_equal(
    unname(vapply(jack_list, function(x) x$estimate, numeric(1))),
    unlist(pop_tajimas_d(test_gt))
  )
  # check that leave-one-out estimates match the subsets by recomputing
  # the jackknife from tajima's d on each subset
  for (i_grp in seq_along(jack_list)) {
    pop_gt <- test_gt %>%
      dplyr::ungroup() %>%
      dplyr::filter(population == names(jack_list)[i_grp])
    taj_full <- pop_tajimas_d(pop_gt)
    taj_loo <- vapply(sort(unique(block_id)), function(i_block) {
      pop_gt %>%
        select_loci(which(block_id != i_block)) %>%
        pop_tajimas_d()
    }, numeric(1))
    m <- as.vector(table(block_id))
    expected <- .jackknife_weighted(taj_full, matrix(taj_loo), matrix(m))
    expect_equal(jack_list[[i_grp]], expected)
  }
})

test_that("jackknife_pop_tajimas_d returns the right format", {
  test_gt <- make_jack_test_gt()
  jack_matrix <- jackknife_pop_tajimas_d(test_gt, block_size = 2)
  jack_tidy <- jackknife_pop_tajimas_d(test_gt, type = "tidy", block_size = 2)
  jack_list <- jackknife_pop_tajimas_d(test_gt, type = "list", block_size = 2)
  expect_equal(names(jack_matrix), c("quantity", "pop1", "pop2", "pop3"))
  expect_equal(
    names(jack_tidy),
    c("group", "estimate", "jack_estimate", "se", "z", "n_blocks")
  )
  expect_equal(names(jack_list), c("pop1", "pop2", "pop3"))
  expect_equal(jack_tidy$se, unname(unlist(jack_matrix[3, -1])))
  expect_equal(jack_list$pop2$se, jack_tidy$se[2])
  # ungrouped data returns a single row
  jack_ungrouped <- jackknife_pop_tajimas_d(dplyr::ungroup(test_gt),
    block_size = 2
  )
  expect_equal(nrow(jack_ungrouped), 1)
  expect_equal(
    jack_ungrouped$estimate,
    pop_tajimas_d(dplyr::ungroup(test_gt))
  )
})

test_that("jackknife_pop_tajimas_d counts pseudohaploid allele copies", {
  pseudohaploid_genotypes <- rbind(
    c(0, 2, NA),
    c(1, 2, 0),
    c(0, NA, 2),
    c(1, 0, 0)
  )
  pseudohaploid_meta <- data.frame(
    id = letters[1:4],
    population = c("popA", "popA", "popB", "popB")
  )
  pseudohaploid_loci <- data.frame(
    name = paste0("rs", 1:3),
    chromosome = 1,
    position = 1:3,
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
  pseudohaploid_gt <- gen_tibble(
    x = pseudohaploid_genotypes,
    loci = pseudohaploid_loci,
    indiv_meta = pseudohaploid_meta,
    quiet = TRUE
  ) %>%
    gt_pseudohaploid(test_n_loci = NULL) %>%
    dplyr::group_by(population)

  jack_taj <- jackknife_pop_tajimas_d(
    pseudohaploid_gt,
    type = "list",
    block_size = 1
  )
  pi_by_group <- loci_pi(pseudohaploid_gt, type = "list")
  ploidy_by_group <- pseudohaploid_gt %>%
    dplyr::summarise(n = sum(indiv_ploidy(genotypes))) %>%
    dplyr::pull(.data$n)

  expected <- Map(tajimas_d_from_pi_vec, pi_by_group, ploidy_by_group)
  expect_equal(
    unname(vapply(jack_taj, function(x) x$estimate, numeric(1))),
    unname(unlist(expected))
  )
})
