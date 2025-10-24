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

    expect_error(gen_tibble(
      vcf_path,
      quiet = TRUE,
      backingfile = tempfile(),
      chunk_size = 2,
      verbose = TRUE,
      convertNA = FALSE,
      parser = "vcfR"
    ), "Unsupported via ...")

    expect_true(all.equal(
      show_genotypes(pop_b_vcf_gt2),
      show_genotypes(pop_b_vcf_gt)
    ))
    expect_true(all.equal(show_loci(pop_b_vcf_gt2), show_loci(pop_b_vcf_gt)))

    # check our cpp parser
    pop_b_vcf_fast_gt <-
      gen_tibble(vcf_path,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "cpp"
      )

    expect_error(
      gen_tibble(vcf_path,
        quiet = TRUE,
        backingfile = tempfile(),
        parser = "cpp",
        verbose = TRUE
      ), "extra parameters can only be used with parser = 'vcfR'"
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
    show_loci(pop_a_bed_gt)[, (names(show_loci(pop_a_bed_gt)) %in% c("big_index", "name", "chromosome", "position", "allele_ref"))], # nolint
    show_loci(pop_a_gt)[, (names(show_loci(pop_a_gt)) %in% c("big_index", "name", "chromosome", "position", "allele_ref"))] # nolint
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
    show_loci(pop_b_bed_gt)[, (names(show_loci(pop_b_bed_gt)) %in% c("big_index", "name", "chromosome", "position", "allele_ref"))], # nolint
    show_loci(pop_b_gt)[, (names(show_loci(pop_b_gt)) %in% c("big_index", "name", "chromosome", "position", "allele_ref"))] # nolint
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
        "extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
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
