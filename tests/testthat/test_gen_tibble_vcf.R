# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

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
    expect_true(file.exists(gt_get_file_names(pop_b_vcf_gt)[1]))
    expect_true(file.exists(gt_get_file_names(pop_b_vcf_gt)[2]))

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

    expect_true(file.exists(gt_get_file_names(pop_b_vcf_fast_gt)[1]))
    expect_true(file.exists(gt_get_file_names(pop_b_vcf_fast_gt)[2]))

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


# test coding names as integers to reduce storage
test_that("vcf tests storing names as integers", {
  vcf_path <- system.file("extdata/pop_a.vcf", package = "tidypopgen")
  pop_a_gt <- gen_tibble(
    vcf_path,
    quiet = TRUE,
    parser = "cpp",
    backingfile = tempfile("pop_a_"),
    names_as_int = TRUE
  )
  expect_true(is.integer(show_loci(pop_a_gt)$name))
  vcf_path <- system.file("extdata/pop_b.vcf", package = "tidypopgen")
  pop_b_gt <- gen_tibble(
    vcf_path,
    quiet = TRUE,
    parser = "cpp",
    backingfile = tempfile("pop_b_"),
    names_as_int = TRUE
  )
  combined_gt <- rbind(pop_a_gt, pop_b_gt,
    use_position = TRUE,
    quiet = TRUE
  )
  expect_true(inherits(combined_gt, "gen_tbl"))
  # and the report keeps the integer names of the loci
  report <- rbind_dry_run(pop_a_gt, pop_b_gt, use_position = TRUE, quiet = TRUE)
  expect_identical(report$ref$name, show_loci(pop_a_gt)$name)
  expect_identical(report$target$name, show_loci(pop_b_gt)$name)

  # integer names can not be used to match loci, so loci are matched by
  # position even with use_position = FALSE (the default), with a message
  # (only once for rbind, which also calls rbind_dry_run)
  # (the report is printed with cat(), so we capture it)
  expect_message(
    capture.output(report_auto <- rbind_dry_run(pop_a_gt, pop_b_gt)),
    "loci are matched by chromosome and position"
  )
  expect_identical(report_auto, report)
  msgs <- character(0)
  capture.output(withCallingHandlers(
    combined_auto <- rbind(pop_a_gt, pop_b_gt,
      backingfile = tempfile("auto_")
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  ))
  expect_identical(
    sum(grepl("loci are matched by chromosome and position", msgs)),
    1L
  )
  expect_identical(show_loci(combined_auto), show_loci(combined_gt))
  expect_identical(
    show_genotypes(combined_auto),
    show_genotypes(combined_gt)
  )
  # no message if quiet
  expect_no_message(
    combined_quiet <- rbind(pop_a_gt, pop_b_gt,
      use_position = FALSE,
      quiet = TRUE,
      backingfile = tempfile("quiet_")
    )
  )
  expect_identical(show_loci(combined_quiet), show_loci(combined_gt))
  expect_no_message(
    rbind_dry_run(pop_a_gt, pop_b_gt, use_position = FALSE, quiet = TRUE)
  )

  # the same if only one of the datasets has integer names
  pop_b_chr_gt <- gen_tibble(
    vcf_path,
    quiet = TRUE,
    parser = "cpp",
    backingfile = tempfile("pop_b_chr_")
  )
  expect_true(is.character(show_loci(pop_b_chr_gt)$name))
  for (ref_target in list(
    list(pop_a_gt, pop_b_chr_gt),
    list(pop_b_chr_gt, pop_a_gt)
  )) {
    expect_identical(
      rbind_dry_run(ref_target[[1]], ref_target[[2]], quiet = TRUE),
      rbind_dry_run(ref_target[[1]], ref_target[[2]],
        use_position = TRUE, quiet = TRUE
      )
    )
    mixed_gt <- rbind(ref_target[[1]], ref_target[[2]],
      quiet = TRUE,
      backingfile = tempfile("mixed_")
    )
    mixed_pos_gt <- rbind(ref_target[[1]], ref_target[[2]],
      use_position = TRUE,
      quiet = TRUE,
      backingfile = tempfile("mixed_pos_")
    )
    expect_identical(show_loci(mixed_gt), show_loci(mixed_pos_gt))
    expect_identical(show_genotypes(mixed_gt), show_genotypes(mixed_pos_gt))
    # the merged loci keep the names (and their type) of the reference
    expect_identical(
      typeof(show_loci(mixed_gt)$name),
      typeof(show_loci(ref_target[[1]])$name)
    )
  }
})

test_that("ped files with a trailing blank line can be read", {
  ped_dir <- tempfile("ped_")
  dir.create(ped_dir)
  ext_dir <- system.file("extdata", package = "tidypopgen")
  file.copy(file.path(ext_dir, "pop_a.map"), ped_dir)
  writeLines(
    c(readLines(file.path(ext_dir, "pop_a.ped")), ""),
    file.path(ped_dir, "pop_a.ped")
  )
  ped_gt <- gen_tibble(file.path(ped_dir, "pop_a.ped"),
    backingfile = tempfile(), quiet = TRUE
  )
  orig_gt <- gen_tibble(file.path(ext_dir, "pop_a.ped"),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(ped_gt$id, orig_gt$id)
  expect_identical(show_genotypes(ped_gt), show_genotypes(orig_gt))
})

test_that("truncated packedancestry geno files are caught", {
  geno_dir <- tempfile("geno_")
  dir.create(geno_dir)
  ext_dir <- system.file("extdata", package = "tidypopgen")
  file.copy(file.path(ext_dir, "pop_a.snp"), geno_dir)
  file.copy(file.path(ext_dir, "pop_a.ind"), geno_dir)
  geno_bytes <- readBin(file.path(ext_dir, "pop_a.geno"), "raw", n = 1e5)
  geno_path <- file.path(geno_dir, "pop_a.geno")
  # remove a few bytes, or most of the file (keeping the header)
  for (n_bytes in c(length(geno_bytes) - 10, 48)) {
    writeBin(geno_bytes[seq_len(n_bytes)], geno_path)
    expect_error(
      gen_tibble(geno_path, backingfile = tempfile(), quiet = TRUE),
      "truncated or corrupted"
    )
  }
})

test_that("file extensions are matched ignoring case", {
  # copy the pop_a files, using upper case extensions for the main files
  upper_dir <- tempfile("upper_")
  dir.create(upper_dir)
  ext_dir <- system.file("extdata", package = "tidypopgen")
  for (ext in c("bim", "fam", "map", "snp", "ind")) {
    file.copy(
      file.path(ext_dir, paste0("pop_a.", ext)),
      file.path(upper_dir, paste0("pop_a.", ext))
    )
  }
  for (ext in c("bed", "ped", "geno")) {
    file.copy(
      file.path(ext_dir, paste0("pop_a.", ext)),
      file.path(upper_dir, paste0("pop_a.", toupper(ext)))
    )
  }
  # a gzipped vcf
  gz_con <- gzfile(file.path(upper_dir, "pop_a.VCF.GZ"), "w")
  writeLines(readLines(file.path(ext_dir, "pop_a.vcf")), gz_con)
  close(gz_con)

  for (ext in c("bed", "ped", "geno", "vcf")) {
    lower_gt <- gen_tibble(file.path(ext_dir, paste0("pop_a.", ext)),
      backingfile = tempfile(), quiet = TRUE
    )
    upper_file <- file.path(
      upper_dir,
      paste0("pop_a.", toupper(ext), ifelse(ext == "vcf", ".GZ", ""))
    )
    upper_gt <- gen_tibble(upper_file, backingfile = tempfile(), quiet = TRUE)
    expect_identical(show_genotypes(upper_gt), show_genotypes(lower_gt))
    expect_identical(show_loci(upper_gt), show_loci(lower_gt))
  }
})

test_that("vcf_gt_dosages_cpp reads genotypes as the R code did", {
  # the R code that vcf_gt_dosages_cpp() replaces, for one genotype
  ref_dosage <- function(g, max_ploidy) {
    x <- strsplit(g, "[/|]")[[1]]
    if (length(x) > max_ploidy) {
      return("error")
    }
    if (length(x) > 0 && !anyNA(x) && !any(x == ".")) {
      dosage <- sum(as.numeric(x))
      if (dosage < (max_ploidy + 1)) as.raw(dosage) else "error"
    } else {
      as.raw(max_ploidy + 1)
    }
  }
  gt <- matrix(
    c("0/1", "1|1", "./.", ".", NA, "", "0|0", "./1", "0/", "1"),
    nrow = 5
  )
  res <- vcf_gt_dosages_cpp(gt, max_ploidy = 2)
  expect_identical(res$error, 0L)
  # individuals are rows
  expect_identical(
    res$dosages,
    matrix(vapply(t(gt), ref_dosage, raw(1), max_ploidy = 2), nrow = 2,
      byrow = FALSE
    )
  )
  # polyploids
  gt <- matrix(c("0/1/1/0", "1/1/1/1", "0/./1/1", "0/0/0/0"), nrow = 2)
  res <- vcf_gt_dosages_cpp(gt, max_ploidy = 4)
  expect_identical(res$dosages, matrix(as.raw(c(2, 5, 4, 0)), nrow = 2))
  # too many alleles, a dosage that is too large, and alleles that are not
  # numbers
  expect_identical(vcf_gt_dosages_cpp(matrix("0/1/1"), 2)$error, 1L)
  expect_identical(vcf_gt_dosages_cpp(matrix("2/1"), 2)$error, 1L)
  expect_identical(vcf_gt_dosages_cpp(matrix("a/1"), 2)$error, 2L)
  expect_identical(vcf_gt_dosages_cpp(matrix("0//1"), 3)$error, 2L)
})

test_that("the md5 sum of the FBM is that of its saved rds for VCF inputs", {
  vcf_path <- system.file("extdata/pop_a.vcf", package = "tidypopgen")
  for (parser in c("cpp", "vcfR")) {
    if (parser == "vcfR") skip_if_not_installed("vcfR")
    test_gt <- gen_tibble(vcf_path,
      parser = parser, quiet = TRUE, backingfile = tempfile()
    )
    rds_file <- gt_get_file_names(test_gt)[1]
    expect_identical(
      attr(test_gt$genotypes, "fbm_md5sum"),
      tools::md5sum(rds_file)
    )
  }
})

test_that("a missing chromosome '.' is NA with both parsers", {
  skip_if_not_installed("vcfR")
  # missing chromosomes, a missing ID (named from chromosome and position)
  # and a missing ALT
  vcf_path <- tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.3",
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb",
      "chrB\t343\trs1\tA\tT\t.\t.\t.\tGT\t0/1\t1/1",
      ".\t65\trs2\tC\tG\t.\t.\t.\tGT\t0/0\t0/1",
      "chrB\t5\t.\tG\tC\t.\t.\t.\tGT\t1/1\t./.",
      ".\t70\t.\tT\t.\t.\t.\t.\tGT\t0/0\t0/0"
    ),
    vcf_path
  )
  cpp_gt <- gen_tibble(vcf_path,
    parser = "cpp", backingfile = tempfile(), quiet = TRUE
  )
  vcfr_gt <- gen_tibble(vcf_path,
    parser = "vcfR", backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(
    as.character(show_loci(cpp_gt)$chromosome),
    c("chrB", NA, "chrB", NA)
  )
  expect_identical(levels(show_loci(cpp_gt)$chromosome), "chrB")
  expect_identical(
    show_loci(cpp_gt)$name,
    c("rs1", "rs2", "chrB_5", "NA_70")
  )
  expect_identical(show_loci(cpp_gt), show_loci(vcfr_gt))
  expect_identical(show_genotypes(cpp_gt), show_genotypes(vcfr_gt))
})
