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
    show_loci(test_gt$genotypes) %>% select(c(-big_index, -chr_int)),
    as_tibble(test_loci)
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
    "can not be a valid allele"
  )
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
    "'x' is not a matrix of integers"
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
    "loci must be one of data.frame or tbl"
  )
})

test_that("gen_tibble required id and population", {
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
    "ind_meta does not include the compulsory column 'id"
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
    "indiv_meta must be one of data.frame, tbl, or list"
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
      "there is a mismatch between the number of loci in the genotype",
      "table x and in the loci table"
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
    "loci does not include the compulsory columns"
  )
})

if (rlang::is_installed("vcfR")) {
  test_that("test gentibble from VCF with missingness issue", {
    ########################
    # PLINK BED files
    ########################
    bed_path <- system.file("extdata/pop_b.bed", package = "tidypopgen")
    pop_b_gt <- gen_tibble(bed_path, quiet = TRUE, backingfile = tempfile())

    ########################
    # PLINK VCF files
    ########################
    vcf_path <- system.file("extdata/pop_b.vcf", package = "tidypopgen")
    pop_b_vcf_gt <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      backingfile = tempfile(),
      parser = "vcfR"
    )
    expect_true(all.equal(
      show_genotypes(pop_b_gt),
      show_genotypes(pop_b_vcf_gt)
    ))
    # reload it in chunks
    pop_b_vcf_gt2 <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      backingfile = tempfile(),
      chunk_size = 2,
      parser = "vcfR"
    )
    expect_true(all.equal(
      show_genotypes(pop_b_vcf_gt2),
      show_genotypes(pop_b_vcf_gt)
    ))
    expect_true(all.equal(show_loci(pop_b_vcf_gt2), show_loci(pop_b_vcf_gt)))
    expect_true(is.integer(show_loci(pop_b_vcf_gt2)$chr_int))

    # check our cpp parser
    pop_b_vcf_fast_gt <-
      gen_tibble(vcf_path,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "cpp"
      )

    # debug
    # if (!identical(show_genotypes(pop_b_gt),
    #                      show_genotypes(pop_b_vcf_fast_gt))) {
    #   stop("genotypes are not equal \n",
    #   print(show_genotypes(pop_b_gt)),"\n",
    #   print(show_genotypes(pop_b_vcf_fast_gt)))
    # }
    #
    # for(i in 1:10){
    #
    #   genotypes_info <- paste0(
    #     "Iteration: ", i, "\n",
    #     "Pop B GT: ", show_genotypes(pop_b_gt), "\n",
    #     "Pop B VCF Fast GT: ", show_genotypes(pop_b_vcf_fast_gt)
    #   )
    #   expect_true(all.equal(show_genotypes(pop_b_gt),
    #                   show_genotypes(pop_b_vcf_fast_gt)), info = genotypes_info) #nolint
    # }
    # check loci table against the vcfR parser
    expect_true(all.equal(
      show_loci(pop_b_vcf_gt),
      show_loci(pop_b_vcf_fast_gt)
    ))
    # reload it in chunks
    pop_b_vcf_fast_gt2 <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      backingfile = tempfile(),
      chunk_size = 2,
      parser = "cpp"
    )

    # # debug #nolint start
    # print(show_genotypes(pop_b_vcf_fast_gt2))
    # print(show_genotypes(pop_b_vcf_fast_gt))
    #
    # for(i in 1:10){
    #
    #   genotypes_info <- paste0(
    #     "Iteration: ", i, "\n",
    #     "Pop B GT: ", show_genotypes(pop_b_vcf_fast_gt), "\n",
    #     "Pop B VCF Fast GT: ", show_genotypes(pop_b_vcf_fast_gt2)
    #   )
    #
    #   expect_true(all.equal(show_genotypes(pop_b_vcf_fast_gt2),
    #                   show_genotypes(pop_b_vcf_fast_gt)), info = genotypes_info)
    # } #nolint end

    expect_true(all.equal(
      show_loci(pop_b_vcf_fast_gt2),
      show_loci(pop_b_vcf_fast_gt)
    ))
    expect_true(is.integer(show_loci(pop_b_vcf_fast_gt2)$chr_int))
  })
}

if (rlang::is_installed("vcfR")) {
  test_that("gen_tibble from files", {
    # test invalid file names
    expect_error(
      gen_tibble("non_existent_file.blah", quiet = TRUE),
      "x should be a valid file path pointing"
    )

    expect_error(
      gen_tibble("non_existent_file.vcf", quiet = TRUE),
      "x should be a valid file path pointing"
    )

    expect_error(
      gen_tibble("", quiet = TRUE),
      "x should not be an empty string"
    )


    ########################
    # PLINK BED files
    ########################
    bed_path <- system.file("extdata/pop_a.bed", package = "tidypopgen")
    pop_a_gt <- gen_tibble(bed_path, quiet = TRUE, backingfile = tempfile())
    # now read the dosages created by plink when saving in raw format
    raw_file_pop_a <-
      read.table(
        system.file("extdata/pop_a.raw", package = "tidypopgen"),
        header = TRUE
      )
    mat <- as.matrix(raw_file_pop_a[, 7:ncol(raw_file_pop_a)])
    mat <- unname(mat)
    expect_true(all.equal(mat, show_genotypes(pop_a_gt)))
    ########################
    # PLINK PED files
    ########################
    ped_path <- system.file("extdata/pop_a.ped", package = "tidypopgen")
    pop_a_ped_gt <- gen_tibble(ped_path, quiet = TRUE, backingfile = tempfile())
    # because ref and alt are defined based on which occurs first in a ped,
    # some alleles will be swapped
    equal_geno <- show_genotypes(pop_a_gt) == show_genotypes(pop_a_ped_gt)
    not_equal <- which(!apply(equal_geno, 2, all))
    # check that the alleles for loci that are mismatched are indeed swapped
    expect_true(all(
      show_loci(pop_a_gt)$allele_alt[not_equal] ==
        show_loci(pop_a_ped_gt)$allele_ref[not_equal]
    ))
    # check that the mismatches are all in the homozygotes
    expect_true(all(
      abs(
        show_genotypes(pop_a_gt)[, not_equal] -
          show_genotypes(pop_a_ped_gt)[, not_equal]
      ) %in%
        c(0, 2)
    ))
    ########################
    # PLINK VCF files
    ########################
    vcf_path <- system.file("extdata/pop_a.vcf", package = "tidypopgen")
    pop_a_vcf_gt <-
      gen_tibble(
        vcf_path,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "vcfR"
      )
    expect_true(all.equal(
      show_genotypes(pop_a_gt),
      show_genotypes(pop_a_vcf_gt)
    ))
    # reload it in chunks
    pop_a_vcf_gt2 <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      backingfile = tempfile(),
      chunk_size = 2,
      parser = "vcfR"
    )
    expect_true(all.equal(
      show_genotypes(pop_a_vcf_gt2),
      show_genotypes(pop_a_vcf_gt)
    ))
    expect_true(all.equal(show_loci(pop_a_vcf_gt2), show_loci(pop_a_vcf_gt)))
    expect_true(is.integer(show_loci(pop_a_vcf_gt)$chr_int))

    # check our cpp parser
    pop_a_vcf_fast_gt <-
      gen_tibble(vcf_path,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "cpp"
      )
    expect_true(all.equal(
      show_genotypes(pop_a_gt),
      show_genotypes(pop_a_vcf_fast_gt)
    ))
    # check loci table against the vcfR parser
    expect_true(all.equal(
      show_loci(pop_a_vcf_gt),
      show_loci(pop_a_vcf_fast_gt)
    ))
    # reload it in chunks
    pop_a_vcf_fast_gt2 <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      backingfile = tempfile(),
      chunk_size = 2,
      parser = "cpp"
    )
    expect_true(all.equal(
      show_genotypes(pop_a_vcf_fast_gt2),
      show_genotypes(pop_a_vcf_fast_gt)
    ))
    expect_true(all.equal(
      show_loci(pop_a_vcf_gt),
      show_loci(pop_a_vcf_fast_gt)
    ))
    expect_true(is.integer(show_loci(pop_a_vcf_fast_gt)$chr_int))
  })
}

if (rlang::is_installed("vcfR")) {
  test_that("gen_tibble from files with missingness", {
    ########################
    # PLINK BED files
    ########################
    bed_path <- system.file("extdata/pop_b.bed", package = "tidypopgen")
    pop_b_gt <- gen_tibble(bed_path, quiet = TRUE, backingfile = tempfile())
    # now read the dosages created by plink when saving in raw format
    raw_file_pop_b <-
      read.table(
        system.file("extdata/pop_b.raw", package = "tidypopgen"),
        header = TRUE
      )
    mat <- as.matrix(raw_file_pop_b[, 7:ncol(raw_file_pop_b)])
    mat <- unname(mat)
    expect_true(all.equal(mat, show_genotypes(pop_b_gt)))
    ########################
    # PLINK PED files
    ########################
    ped_path <- system.file("extdata/pop_b.ped", package = "tidypopgen")
    pop_b_ped_gt <- gen_tibble(ped_path, quiet = TRUE, backingfile = tempfile())
    # because ref and alt are defined based on which occurs first in a ped,
    # some alleles will be swapped
    equal_geno <- show_genotypes(pop_b_gt) == show_genotypes(pop_b_ped_gt)
    not_equal <- which(!apply(equal_geno, 2, all))
    # check that the alleles for loci that are mismatched are indeed swapped
    expect_true(all(
      show_loci(pop_b_gt)$allele_alt[not_equal] ==
        show_loci(pop_b_ped_gt)$allele_ref[not_equal]
    ))
    # check that the mismatches are all in the homozygotes
    expect_true(all(
      abs(
        show_genotypes(pop_b_gt)[, not_equal] -
          show_genotypes(pop_b_ped_gt)[, not_equal]
      ) %in%
        c(0, 2)
    ))
    ########################
    # PLINK VCF files
    ########################
    vcf_path <- system.file("extdata/pop_b.vcf", package = "tidypopgen")
    pop_b_vcf_gt <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      backingfile = tempfile(),
      parser = "vcfR"
    )
    expect_true(all.equal(
      show_genotypes(pop_b_gt),
      show_genotypes(pop_b_vcf_gt)
    ))
    # reload it in chunks
    pop_b_vcf_gt2 <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      backingfile = tempfile(),
      chunk_size = 2,
      parser = "vcfR"
    )
    expect_true(all.equal(
      show_genotypes(pop_b_vcf_gt2),
      show_genotypes(pop_b_vcf_gt)
    ))
    expect_true(all.equal(show_loci(pop_b_vcf_gt2), show_loci(pop_b_vcf_gt)))
    expect_true(is.integer(show_loci(pop_b_vcf_gt2)$chr_int))
  })
}

test_that("gentibble with packedancestry", {
  geno_path <- system.file("extdata/pop_a.geno", package = "tidypopgen")
  pop_a_gt <-
    gen_tibble(
      geno_path,
      quiet = TRUE,
      backingfile = tempfile(),
      valid_alleles = c("A", "G", "C", "T")
    )
  # dosages in packedancestry files are the opposite to plink files
  # our new packedancestry reader will flip dosages and alleles
  # therefore a gt from packedancestry should now match a gt from .bed

  # The order our files were created in matters here, because plink
  # swaps alleles based on frequency
  # pop_a.bed was created first and used to make pop_a.raw and pop_a.ped files
  # pop_a.ped was then used to create pop_a.geno

  ########################
  # Compare to PED files
  ########################
  ped_path <- system.file("extdata/pop_a.ped", package = "tidypopgen")
  pop_a_ped_gt <- gen_tibble(ped_path, quiet = TRUE, backingfile = tempfile())
  # because ref and alt are defined based on which occurs first in a ped,
  # some alleles will be swapped
  equal_geno <- show_genotypes(pop_a_gt) == show_genotypes(pop_a_ped_gt)
  not_equal <- which(!apply(equal_geno, 2, all))
  # check that the alleles for loci that are mismatched are indeed swapped
  expect_true(all(
    show_loci(pop_a_gt)$allele_alt[not_equal] ==
      show_loci(pop_a_ped_gt)$allele_ref[not_equal]
  ))
  # check that the mismatches are all in the homozygotes
  expect_true(all(
    abs(
      show_genotypes(pop_a_gt)[, not_equal] -
        show_genotypes(pop_a_ped_gt)[, not_equal]
    ) %in%
      c(0, 2)
  ))

  ########################
  # Compare to .raw
  ########################
  raw_file_pop_a <-
    read.table(
      system.file("extdata/pop_a.raw", package = "tidypopgen"),
      header = TRUE
    )
  mat <- as.matrix(raw_file_pop_a[, 7:ncol(raw_file_pop_a)])
  mat <- unname(mat)
  expect_true(all.equal(mat, show_genotypes(pop_a_gt)))

  ########################
  # Compare to BED
  ########################
  bed_path <- system.file("extdata/pop_a.bed", package = "tidypopgen")
  pop_a_bed_gt <- gen_tibble(bed_path, quiet = TRUE, backingfile = tempfile())
  # we expect the genotypes to be exactly the same
  expect_true(all.equal(show_genotypes(pop_a_bed_gt), show_genotypes(pop_a_gt)))
  # we expect the loci info to be the same apart from genetic.dist, which is not
  # recorded, and allele_alt, because .geno was written from .ped, and when data
  # are loaded into a gen_tibble from .ped, allele_alt is NA when all
  # individuals are homozygous ref. This is the case for snps rs9697457,
  # rs2862633, rs28569024.
  expect_true(all.equal(
    show_loci(pop_a_bed_gt)[, (names(show_loci(pop_a_bed_gt)) %in% c("big_index", "name", "chromosome", "position", "allele_ref", "chr_int"))], # nolint
    show_loci(pop_a_gt)[, (names(show_loci(pop_a_gt)) %in% c("big_index", "name", "chromosome", "position", "allele_ref", "chr_int"))] # nolint
  ))
})

test_that("gentibble with packedancestry and missingness", {
  geno_path <- system.file("extdata/pop_b.geno", package = "tidypopgen")
  pop_b_gt <-
    gen_tibble(
      geno_path,
      quiet = TRUE,
      backingfile = tempfile(),
      valid_alleles = c("A", "G", "C", "T")
    )
  # dosages in packedancestry files are the opposite to plink files
  # our new packedancestry reader will flip dosages and alleles
  # therefore a gt from packedancestry should now match a gt from .bed

  # The order our files were created in matters here, because plink
  # swaps alleles based on frequency
  # pop_a.bed was created first and used to make pop_a.raw and pop_a.ped files
  # pop_a.ped was then used to create pop_a.geno

  ########################
  # Compare to PED files
  ########################
  ped_path <- system.file("extdata/pop_b.ped", package = "tidypopgen")
  pop_b_ped_gt <- gen_tibble(ped_path, quiet = TRUE, backingfile = tempfile())
  # because ref and alt are defined based on which occurs first in a ped,
  # some alleles will be swapped
  equal_geno <- show_genotypes(pop_b_gt) == show_genotypes(pop_b_ped_gt)
  not_equal <- which(!apply(equal_geno, 2, all))
  # check that the alleles for loci that are mismatched are indeed swapped
  expect_true(all(
    show_loci(pop_b_gt)$allele_alt[not_equal] ==
      show_loci(pop_b_ped_gt)$allele_ref[not_equal]
  ))
  # check that the mismatches are all in the homozygotes
  expect_true(all(
    abs(
      show_genotypes(pop_b_gt)[, not_equal] -
        show_genotypes(pop_b_ped_gt)[, not_equal]
    ) %in%
      c(0, 2)
  ))

  ########################
  # Compare to .raw
  ########################
  raw_file_pop_b <-
    read.table(
      system.file("extdata/pop_b.raw", package = "tidypopgen"),
      header = TRUE
    )
  mat <- as.matrix(raw_file_pop_b[, 7:ncol(raw_file_pop_b)])
  mat <- unname(mat)
  expect_true(all.equal(mat, show_genotypes(pop_b_gt)))

  ########################
  # Compare to BED
  ########################
  bed_path <- system.file("extdata/pop_b.bed", package = "tidypopgen")
  pop_b_bed_gt <- gen_tibble(bed_path, quiet = TRUE, backingfile = tempfile())
  # we expect the genotypes to be exactly the same
  expect_true(all.equal(show_genotypes(pop_b_bed_gt), show_genotypes(pop_b_gt)))
  # we expect the loci info to be the same apart from genetic.dist, which is not
  # recorded, and allele_alt, because .geno was written from .ped, and when data
  # are loaded into a gen_tibble from .ped, allele_alt is NA when all
  # individuals are homozygous ref. This is the case for snps rs9697457,
  # rs2862633, rs28569024.
  expect_true(all.equal(
    show_loci(pop_b_bed_gt)[, (names(show_loci(pop_b_bed_gt)) %in% c("big_index", "name", "chromosome", "position", "allele_ref", "chr_int"))], # nolint
    show_loci(pop_b_gt)[, (names(show_loci(pop_b_gt)) %in% c("big_index", "name", "chromosome", "position", "allele_ref", "chr_int"))] # nolint
  ))
})

test_that("check summary stats for gen_tibbles read in different ways", {
  geno_path <- system.file("extdata/pop_a.geno", package = "tidypopgen")
  pop_a_gt_geno <-
    gen_tibble(
      geno_path,
      quiet = TRUE,
      backingfile = tempfile(),
      valid_alleles = c("A", "G", "C", "T")
    )

  bed_path <- system.file("extdata/pop_a.bed", package = "tidypopgen")
  pop_a_gt_bed <-
    gen_tibble(
      bed_path,
      quiet = TRUE,
      backingfile = tempfile(),
      valid_alleles = c("A", "G", "C", "T")
    )

  ped_path <- system.file("extdata/pop_a.ped", package = "tidypopgen")
  pop_a_gt_ped <-
    gen_tibble(
      ped_path,
      quiet = TRUE,
      backingfile = tempfile(),
      valid_alleles = c("A", "G", "C", "T")
    )

  # MAF
  bed_maf <- pop_a_gt_bed %>% loci_maf()
  ped_maf <- pop_a_gt_ped %>% loci_maf()
  geno_maf <- pop_a_gt_geno %>% loci_maf()
  expect_equal(bed_maf, ped_maf)
  expect_equal(bed_maf, geno_maf)
  expect_equal(ped_maf, geno_maf)

  # Missing loci
  bed_miss <- pop_a_gt_bed %>% loci_missingness()
  ped_miss <- pop_a_gt_ped %>% loci_missingness()
  geno_miss <- pop_a_gt_geno %>% loci_missingness()

  expect_equal(bed_miss, ped_miss)
  expect_equal(bed_miss, geno_miss)
  expect_equal(ped_miss, geno_miss)

  # Missing individuals
  bed_miss <- pop_a_gt_bed %>% indiv_missingness()
  ped_miss <- pop_a_gt_ped %>% indiv_missingness()
  geno_miss <- pop_a_gt_geno %>% indiv_missingness()

  expect_equal(bed_miss, ped_miss)
  expect_equal(bed_miss, geno_miss)
  expect_equal(ped_miss, geno_miss)
})

if (rlang::is_installed("vcfR")) {
  test_that("chr_int is always an integer", {
    # matrix method
    test_gt <- gen_tibble(
      x = test_genotypes,
      loci = test_loci,
      indiv_meta = test_indiv_meta,
      quiet = TRUE
    )
    expect_true(is.integer(show_loci(test_gt)$chr_int))

    test_loci_fac <- data.frame(
      name = paste0("rs", 1:6),
      chromosome = as.factor(paste0("chr", c(1, 1, 1, 1, 2, 2))),
      position = as.integer(c(3, 5, 65, 343, 23, 456)),
      genetic_dist = as.double(rep(0, 6)),
      allele_ref = c("A", "T", "C", "G", "C", "T"),
      allele_alt = c("T", "C", NA, "C", "G", "A")
    )
    test_gt <- gen_tibble(
      x = test_genotypes,
      loci = test_loci_fac,
      indiv_meta = test_indiv_meta,
      quiet = TRUE
    )
    expect_true(is.integer(show_loci(test_gt)$chr_int))

    # character methods
    bed_path <- system.file("extdata/pop_a.bed", package = "tidypopgen")
    pop_a_gt <- gen_tibble(bed_path, quiet = TRUE, backingfile = tempfile())
    expect_true(is.integer(show_loci(pop_a_gt)$chr_int))

    ped_path <- system.file("extdata/pop_a.ped", package = "tidypopgen")
    pop_a_ped_gt <- gen_tibble(ped_path, quiet = TRUE, backingfile = tempfile())
    expect_true(is.integer(show_loci(pop_a_ped_gt)$chr_int))

    vcf_path <- system.file("extdata/pop_a.vcf", package = "tidypopgen")
    pop_a_vcf_gt <-
      gen_tibble(
        vcf_path,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "vcfR"
      )
    expect_true(is.integer(show_loci(pop_a_vcf_gt)$chr_int))
    pop_a_vcf_fast_gt <-
      gen_tibble(vcf_path,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "cpp"
      )
    expect_true(is.integer(show_loci(pop_a_vcf_fast_gt)$chr_int))

    geno_path <- system.file("extdata/pop_a.geno", package = "tidypopgen")
    pop_a_gt <-
      gen_tibble(
        geno_path,
        quiet = TRUE,
        backingfile = tempfile(),
        valid_alleles = c("A", "G", "C", "T")
      )
    expect_true(is.integer(show_loci(pop_a_gt)$chr_int))
    geno_path <- system.file("extdata/pop_b.geno", package = "tidypopgen")
    pop_b_gt <-
      gen_tibble(
        geno_path,
        quiet = TRUE,
        backingfile = tempfile(),
        valid_alleles = c("A", "G", "C", "T")
      )
    expect_true(is.integer(show_loci(pop_b_gt)$chr_int))
  })
}

test_that("gt without population is valid", {
  test_indiv_meta <- data.frame(id = c("a", "b", "c"))
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
    quiet = TRUE,
    backingfile = tempfile()
  )

  # still inherits gen_tbl
  stopifnot_gen_tibble(test_gt)
  expect_true(inherits(test_gt, "gen_tbl"))

  # can save an load as usual
  file_names <- gt_save(test_gt, file = tempfile(), quiet = TRUE)
  gt_load(file_names[1])

  # will fail in functions that require grouping
  expect_error(pop_fst(test_gt), ".x should be a grouped gen_tibble")
  # functions that require grouping work after adding groups
  test_gt$groups <- c("A", "A", "B")
  test_gt <- test_gt %>% group_by(groups)
  expect_equal(unname(pop_fst(test_gt)), c(-0, NaN))

  # gen_tbl from vcf does not have population automatically
  vcf_path <- system.file("extdata/pop_b.vcf", package = "tidypopgen")
  pop_b_vcf_gt <- gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile())
  expect_true(!("population" %in% names(pop_b_vcf_gt)))
})

if (rlang::is_installed("vcfR")) {
  test_that("additional vcf tests with larger file", {
    vcf_path <-
      system.file(
        "/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
        package = "tidypopgen"
      )
    anole_gt <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      parser = "cpp",
      backingfile = tempfile("anolis_")
    )
    anole_gt_vcfr <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      parser = "vcfR",
      backingfile = tempfile("anolis_")
    )
    expect_true(all.equal(
      show_genotypes(anole_gt),
      show_genotypes(anole_gt_vcfr)
    ))
    anole_gt2 <- gen_tibble(
      vcf_path,
      quiet = TRUE,
      parser = "cpp",
      backingfile = tempfile("anolis_"),
      chunk_size = 1000,
      n_cores = 2
    )
    expect_true(all.equal(
      show_genotypes(anole_gt2),
      show_genotypes(anole_gt_vcfr)
    ))
  })
}

if (rlang::is_installed("vcfR")) {
  test_that("vcf's with haploid markers", {
    # read vcf with haploid markers first
    vcf_path_haploid <- system.file(
      "extdata/haploid_first_pop_a.vcf",
      package = "tidypopgen"
    )
    # the cpp parser catches the problem
    expect_error(
      pop_a_vcf_gt_hap_cpp <-
        gen_tibble(
          vcf_path_haploid,
          quiet = TRUE,
          backingfile = tempfile(),
          parser = "cpp"
        ),
      "a genotype"
    )
    # vcfR catches the problem
    expect_error(
      pop_a_vcf_gt_hap_vcfr <-
        gen_tibble(
          vcf_path_haploid,
          quiet = TRUE,
          backingfile = tempfile(),
          parser = "vcfR"
        ),
      "a genotype"
    )
    # read vcf with haploid in the middle
    vcf_path_haploid_middle <-
      system.file("extdata/haploid_middle_pop_a.vcf", package = "tidypopgen")
    pop_a_vcf_gt_hapmid_cpp <-
      gen_tibble(
        vcf_path_haploid_middle,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "cpp"
      )
    pop_a_vcf_gt_hapmid_vcfr <-
      gen_tibble(
        vcf_path_haploid_middle,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "vcfR"
      )

    # vcfr reads correctly
    expect_equal(
      show_genotypes(pop_a_vcf_gt_hapmid_vcfr)[
        ,
        show_loci(pop_a_vcf_gt_hapmid_vcfr)$chromosome == 23
      ],
      c(1, 0, 0, 0, 0)
    )
    # cpp reads correctly
    expect_equal(
      show_genotypes(pop_a_vcf_gt_hapmid_cpp)[
        ,
        show_loci(pop_a_vcf_gt_hapmid_cpp)$chromosome == 23
      ],
      c(1, 0, 0, 0, 0)
    )
  })
}

if (rlang::is_installed("vcfR")) {
  test_that("chr_int is correct", {
    # unit tests for the casting function
    chromosome_names <- c("1", "2", NA, "4")
    expect_true(identical(
      c(1L, 2L, NA, 4L),
      cast_chromosome_to_int(chromosome_names)
    ))
    chromosome_names <- c("chr1", "chr2", NA, "chr4")
    expect_true(identical(
      c(1L, 2L, NA, 4L),
      cast_chromosome_to_int(chromosome_names)
    ))
    chromosome_names <- c("a", "b", NA, "c")
    expect_true(identical(
      c(1L, 2L, NA, 3L),
      cast_chromosome_to_int(chromosome_names)
    ))

    # a real life example

    # read bed
    bed_path <- system.file("extdata/pop_b.bed", package = "tidypopgen")
    pop_b_gt <- gen_tibble(bed_path, quiet = TRUE, backingfile = tempfile())

    # chr_int is correct for bed files
    expect_equal(
      as.character(show_loci(pop_b_gt)$chr_int),
      show_loci(pop_b_gt)$chromosome
    )

    # read ped
    ped_path <- system.file("extdata/pop_b.ped", package = "tidypopgen")
    pop_b_ped_gt <- gen_tibble(ped_path, quiet = TRUE, backingfile = tempfile())

    # chr_int is correct for ped files
    expect_equal(
      as.character(show_loci(pop_b_ped_gt)$chr_int),
      show_loci(pop_b_ped_gt)$chromosome
    )

    # read vcf
    vcf_path <- system.file("extdata/pop_a.vcf", package = "tidypopgen")
    pop_a_vcfr_gt <-
      gen_tibble(
        vcf_path,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "vcfR"
      )
    pop_a_cpp_gt <-
      gen_tibble(vcf_path,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "cpp"
      )

    # chr_int is correct for vcf
    expect_equal(
      show_loci(pop_a_vcfr_gt)$chr_int,
      as.integer(show_loci(pop_a_vcfr_gt)$chromosome)
    )
    expect_equal(
      show_loci(pop_a_cpp_gt)$chr_int,
      as.integer(show_loci(pop_a_cpp_gt)$chromosome)
    )
  })
}

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
