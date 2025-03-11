raw_path_pop_a <- system.file("extdata/pop_a.bed", package = "tidypopgen")
bigsnp_path_a <- bigsnpr::snp_readBed(
  raw_path_pop_a,
  backingfile = tempfile("test_a_")
)
pop_a_snp <- bigsnpr::snp_attach(bigsnp_path_a)

geno_start <- pop_a_snp$genotypes[]

# quick function to swap loci
swap_locus <- function(x) {
  x_new <- x
  x_new[x == 2] <- 0
  x_new[x == 0] <- 2
  x_new
}

test_that("subset_bigSNP correctly subsets the data", {
  indiv_indices <- c(2, 1, 4)
  loci_indices <- c(1:5, 6, 11, 10)
  swap_indices <- c(2, 5)
  new_snp <- subset_bigSNP(
    pop_a_snp,
    indiv_indices = indiv_indices,
    loci_indices = loci_indices,
    swap_indices = swap_indices
  )
  # we first swap the old loci
  geno_start[, swap_indices] <- swap_locus(geno_start[, swap_indices])
  expect_identical(new_snp$genotypes[], geno_start[indiv_indices, loci_indices])
  # check the swapped locus map
  swapped_locus1 <- pop_a_snp$map$marker.ID[swap_indices[1]]
  expect_true(
    pop_a_snp$map %>%
      dplyr::filter(marker.ID == swapped_locus1) %>%
      pull(allele1) ==
      new_snp$map %>%
        dplyr::filter(marker.ID == swapped_locus1) %>%
        pull(allele2)
  )

  ##############################################################################
  ## check missing data
  # add some missing data
  pop_a_snp$genotypes[1, 5] <- 3
  # check that this was successful
  expect_true(is.na(pop_a_snp$genotypes[1, 5]))
  # now do the subset with swap
  new_snp2 <- subset_bigSNP(
    pop_a_snp,
    indiv_indices = indiv_indices,
    loci_indices = loci_indices,
    swap_indices = swap_indices
  )
  # check that we preserved the NA
  expect_true(is.na(new_snp2$genotypes[2, 5]))
  # this should be identical to the original subset, except for the na value
  new_geno_copy <- new_snp$genotypes[]
  new_geno_copy[2, 5] <- NA
  expect_true(all.equal(new_snp2$genotypes[], new_geno_copy))
  # now impute
  pop_a_snp$genotypes <- bigsnpr::snp_fastImputeSimple(pop_a_snp$genotypes)
  # check that it worked
  expect_false(is.na(pop_a_snp$genotypes[1, 5]))
  # now subset like before
  new_snp3 <- subset_bigSNP(
    pop_a_snp,
    indiv_indices = indiv_indices,
    loci_indices = loci_indices,
    swap_indices = swap_indices
  )
  # check that we don't have an NA
  expect_false(is.na(new_snp3$genotypes[2, 5]))
  # and that we flipped it correctly
  # (it should be identical to the original subset)
  expect_true(all.equal(new_snp3$genotypes[, 5], new_snp$genotypes[, 5]))
  # but if we change the code, it is an imputed genotype
  new_snp3$genotypes$code256 <- bigsnpr::CODE_012
  expect_true(is.na(new_snp3$genotypes[2, 5]))
})
