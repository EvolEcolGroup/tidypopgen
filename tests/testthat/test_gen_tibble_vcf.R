# delete the backing files of this file when it ends (see helper_cleanup.R)
local_tempdir_cleanup()

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

test_that("order_as_levels sorts as the levels of as.factor()", {
  # (unique strings, as the chromosome names given by the cpp parser)
  examples <- list(
    paste0("SNP_higher_path_", c(10, 2, 1, 100)),
    c(paste0("chr", c(10, 2, 1)), "chrX", "chrUn_KI270302v1"),
    # mixed case and punctuation, which C order and collation can sort
    # differently
    c("chr_1", "chr1", "Chr2", "chr10", "a", "B", "b"),
    "1",
    character(0)
  )
  check_levels <- function() {
    for (x in examples) {
      expect_identical(x[order_as_levels(x)], levels(as.factor(x)))
    }
  }
  # in the collation of the tests (C), and in an English locale if available
  # (checked by setting it, and then restoring the collation)
  check_levels()
  old_collate <- Sys.getlocale("LC_COLLATE")
  has_english <- suppressWarnings(
    Sys.setlocale("LC_COLLATE", "en_US.UTF-8")
  ) != ""
  Sys.setlocale("LC_COLLATE", old_collate)
  if (has_english) {
    withr::local_collate("en_US.UTF-8")
    check_levels()
  }
})

test_that("malformed variant lines give clear errors with the cpp parser", {
  header <- c(
    "##fileformat=VCFv4.3",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb",
    "chr1\t5\trs1\tA\tT\t.\t.\t.\tGT\t0/1\t1/1"
  )
  # a line with fewer than 5 fields
  short_vcf <- tempfile(fileext = ".vcf")
  writeLines(c(header, "chr1\t9"), short_vcf)
  expect_error(
    gen_tibble(short_vcf,
      parser = "cpp", backingfile = tempfile(), quiet = TRUE
    ),
    "fewer than 5 fields"
  )
  # a position that is not a number
  bad_pos_vcf <- tempfile(fileext = ".vcf")
  writeLines(
    c(header, "chr1\t.\trs2\tC\tG\t.\t.\t.\tGT\t0/0\t0/1"),
    bad_pos_vcf
  )
  expect_error(
    gen_tibble(bad_pos_vcf,
      parser = "cpp", backingfile = tempfile(), quiet = TRUE
    ),
    "not a valid integer"
  )
})

test_that("the cpp parser codes missing genotypes once ploidies are known", {
  # b is missing in the first variants, and has a larger ploidy (4) than the
  # others, so its ploidy (and the value for missing genotypes, max_ploidy + 1)
  # is only known after some genotypes have been read
  vcf_path <- tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.3",
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc",
      "chr1\t5\trs1\tA\tT\t.\t.\t.\tGT\t0/1\t.\t1/1",
      "chr1\t7\trs2\tA\tT\t.\t.\t.\tGT\t./.\t.\t0/1",
      "chr1\t9\trs3\tC\tG\t.\t.\t.\tGT\t0/0\t0/1/1/1\t./.",
      "chr1\t12\trs4\tG\tA\t.\t.\t.\tGT\t1/1\t1/1/0/0\t0/1"
    ),
    vcf_path
  )
  test_gt <- gen_tibble(vcf_path,
    parser = "cpp", backingfile = tempfile(), quiet = TRUE
  )
  expect_equal(indiv_ploidy(test_gt), c(2, 4, 2))
  expect_equal(
    show_genotypes(test_gt),
    rbind(c(1, NA, 0, 2), c(NA, NA, 3, 2), c(2, 1, NA, 1)),
    ignore_attr = TRUE
  )
  # the backing file stores missing genotypes as max_ploidy + 1 (5)
  bk_bytes <- readBin(gt_get_file_names(test_gt)[2], what = "raw", n = 12)
  expect_identical(
    as.integer(bk_bytes),
    c(1L, 5L, 2L, 5L, 5L, 1L, 0L, 3L, 5L, 2L, 2L, 1L)
  )
})

test_that("the cpp parser leaves no backing files after an error", {
  # a genotype with more alleles than estimated from the first variant
  backingfile <- tempfile()
  expect_error(
    gen_tibble(
      system.file("extdata/haploid_first_pop_a.vcf", package = "tidypopgen"),
      parser = "cpp", backingfile = backingfile, quiet = TRUE
    ),
    "a genotype has more than max_ploidy alleles"
  )
  expect_length(
    list.files(dirname(backingfile), pattern = basename(backingfile)),
    0
  )
  # a VCF with no biallelic SNPs (an indel and a multiallelic SNP)
  vcf_path <- tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.3",
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb",
      "chr1\t5\trs1\tAT\tA\t.\t.\t.\tGT\t0/1\t0/0",
      "chr1\t9\trs2\tC\tG,T\t.\t.\t.\tGT\t0/0\t0/2"
    ),
    vcf_path
  )
  backingfile <- tempfile()
  expect_error(
    gen_tibble(vcf_path,
      parser = "cpp", backingfile = backingfile, quiet = TRUE
    ),
    "no biallelic SNPs"
  )
  expect_length(
    list.files(dirname(backingfile), pattern = basename(backingfile)),
    0
  )
  # a VCF with more genotypes than samples in the header (an extra tab)
  writeLines(
    c(
      "##fileformat=VCFv4.3",
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb",
      "chr1\t5\trs1\tA\tT\t.\t.\t.\tGT\t0/1\t\t1/1"
    ),
    vcf_path
  )
  backingfile <- tempfile()
  expect_error(
    gen_tibble(vcf_path,
      parser = "cpp", backingfile = backingfile, quiet = TRUE
    ),
    "has 3 genotypes, but the header line \\(#CHROM\\) has 2 samples"
  )
  expect_length(
    list.files(dirname(backingfile), pattern = basename(backingfile)),
    0
  )
})

test_that("the cpp parser reads bgzipped VCFs in parallel", {
  # this VCF is compressed with bgzip (BGZF, made of blocks of at most 64 KB,
  # which are decompressed in parallel)
  bgzf_path <- system.file(
    "extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
    package = "tidypopgen"
  )
  bgzf_gt <- gen_tibble(bgzf_path,
    parser = "cpp", n_cores = 2, backingfile = tempfile(), quiet = TRUE
  )
  # the same VCF compressed with gzip (as a single stream, which is
  # decompressed sequentially)
  gzip_path <- tempfile(fileext = ".vcf.gz")
  gz_con <- gzfile(gzip_path, "w")
  writeLines(readLines(bgzf_path), gz_con)
  close(gz_con)
  gzip_gt <- gen_tibble(gzip_path,
    parser = "cpp", n_cores = 1, backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(bgzf_gt$id, gzip_gt$id)
  expect_identical(show_loci(bgzf_gt), show_loci(gzip_gt))
  expect_identical(show_genotypes(bgzf_gt), show_genotypes(gzip_gt))
  # a corrupted block gives an error (its CRC32 does not match)
  bgzf_bytes <- readBin(bgzf_path, what = "raw", n = file.size(bgzf_path))
  middle <- length(bgzf_bytes) %/% 2
  bgzf_bytes[middle] <- xor(bgzf_bytes[middle], as.raw(255))
  corrupt_path <- tempfile(fileext = ".vcf.gz")
  writeBin(bgzf_bytes, corrupt_path)
  expect_error(
    gen_tibble(corrupt_path,
      parser = "cpp", n_cores = 2, backingfile = tempfile(), quiet = TRUE
    ),
    "error decompressing the BGZF file"
  )
  # the first block of the file (whose size, minus 1, is in bytes 17-18 of
  # its header), truncated or followed by bytes that are not a BGZF block
  block_size <- sum(as.integer(bgzf_bytes[17:18]) * c(1L, 256L)) + 1L
  truncated_path <- tempfile(fileext = ".vcf.gz")
  writeBin(bgzf_bytes[seq_len(block_size - 10)], truncated_path)
  expect_error(
    gen_tibble(truncated_path, backingfile = tempfile(), quiet = TRUE),
    "the BGZF file is truncated"
  )
  junk_path <- tempfile(fileext = ".vcf.gz")
  writeBin(c(bgzf_bytes[seq_len(block_size)], charToRaw("junk")), junk_path)
  expect_error(
    gen_tibble(junk_path, backingfile = tempfile(), quiet = TRUE),
    "not a valid BGZF file"
  )
  # a block whose uncompressed size (its last 4 bytes) is larger than a BGZF
  # block can be
  first_block <- bgzf_bytes[seq_len(block_size)]
  first_block[block_size - 0:3] <- as.raw(255)
  big_path <- tempfile(fileext = ".vcf.gz")
  writeBin(first_block, big_path)
  expect_error(
    gen_tibble(big_path, backingfile = tempfile(), quiet = TRUE),
    "error decompressing the BGZF file"
  )
})

test_that("the cpp parser rejects invalid genotypes and positions", {
  loci <- data.frame(
    chrom = "1", pos = c(100, 200), id = c("s1", "s2"),
    ref = c("A", "G"), alt = c("C", "T")
  )
  # a VCF with the genotypes geno at its second locus (with FORMAT `format`,
  # and position `pos`)
  bad_vcf <- function(geno, format = "GT", pos = "200") {
    vcf_path <- write_merge_vcf(c("a", "b"), loci,
      genotypes = rbind(c("0/1", "0/0"), geno)
    )
    vcf_lines <- readLines(vcf_path)
    vcf_lines[5] <- sub("\tGT\t", paste0("\t", format, "\t"), vcf_lines[5])
    vcf_lines[5] <- sub("\t200\t", paste0("\t", pos, "\t"), vcf_lines[5])
    writeLines(vcf_lines, vcf_path)
    vcf_path
  }
  # reading it gives an error, and leaves no backing file
  expect_vcf_error <- function(vcf_path, regexp) {
    backingfile <- tempfile()
    expect_error(
      gen_tibble(vcf_path, backingfile = backingfile, quiet = TRUE),
      regexp
    )
    expect_length(
      list.files(dirname(backingfile), pattern = basename(backingfile)),
      0
    )
  }
  # these genotypes are fine (values after the genotype are ignored)
  good_gt <- gen_tibble(bad_vcf(c("1/1", "0|1:20")),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_equal(show_genotypes(good_gt), cbind(c(1, 0), c(2, 1)))
  # alleles that a biallelic SNP can not have, an allele with two digits, or
  # an empty genotype
  for (geno in list(c("2/1", "0/0"), c("a/1", "0/0"), c("0/12", "0/0"),
                    c("1/1", ":20"), c("1/1", ""))) {
    expect_vcf_error(bad_vcf(geno), "is not valid: the alleles of a")
  }
  # the genotype is not the first value of each sample
  expect_vcf_error(
    bad_vcf(c("20:1/1", "20:0/0"), format = "DP:GT"),
    "does not start with GT"
  )
  # a position that is not an integer
  expect_vcf_error(bad_vcf(c("1/1", "0/0"), pos = "200x"), "valid integer")
  # an invalid number of cores
  vcf_path <- write_merge_vcf(c("a", "b"), loci,
    genotypes = rbind(c("0/1", "0/0"), c("1/1", "0/1"))
  )
  for (n_cores in list(0, NA, Inf, "2")) {
    expect_error(
      gen_tibble(vcf_path, n_cores = n_cores, backingfile = tempfile()),
      "'n_cores' should be a positive integer"
    )
  }
  expect_error(
    tidypopgen:::vcf_to_bk(vcf_path, bk_file = "", n_threads = 0,
      read_genotypes = FALSE
    ),
    "'n_threads' should be a positive integer"
  )
})

test_that("the cpp parser handles unusual but valid VCF files", {
  loci <- data.frame(
    chrom = "1", pos = c(100, 200), id = c("s1", "s2"),
    ref = c("A", "G"), alt = c("C", "T")
  )
  vcf_path <- write_merge_vcf(c("a", "b"), loci,
    genotypes = rbind(c("0/1", "0/0"), c("1/1", "0/1"))
  )
  vcf_lines <- readLines(vcf_path)
  # a file without a newline at the end of its last line
  no_newline_path <- tempfile(fileext = ".vcf")
  cat(paste(vcf_lines, collapse = "\n"), file = no_newline_path)
  no_newline_gt <- gen_tibble(no_newline_path,
    backingfile = tempfile(), quiet = TRUE
  )
  expect_equal(show_genotypes(no_newline_gt), rbind(c(1, 2), c(0, 1)))
  # a sample with only missing genotypes is diploid
  missing_gt <- gen_tibble(
    write_merge_vcf("a", loci, genotypes = rbind(".", ".")),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(indiv_ploidy(missing_gt), 2L)
  expect_equal(show_genotypes(missing_gt), rbind(c(NA_real_, NA_real_)))
  # a variant line (not the first) with too few genotypes is an error
  n_lines <- length(vcf_lines)
  vcf_lines[n_lines] <- sub("\t[^\t]*$", "", vcf_lines[n_lines])
  short_path <- tempfile(fileext = ".vcf")
  writeLines(vcf_lines, short_path)
  expect_error(
    gen_tibble(short_path, backingfile = tempfile(), quiet = TRUE),
    "fewer genotypes than individuals"
  )
})

test_that("vcf_to_fbm_cpp puts the backing file next to the VCF", {
  vcf_dir <- tempfile()
  dir.create(vcf_dir)
  vcf_path <- file.path(vcf_dir, "small.vcf")
  write_merge_vcf("a",
    data.frame(chrom = "1", pos = 100, id = "s1", ref = "A", alt = "C"),
    genotypes = rbind("0/1"), file = vcf_path
  )
  small_gt <- vcf_to_fbm_cpp(vcf_path, quiet = TRUE)
  expect_identical(
    basename(attr(small_gt$genotypes, "fbm")$backingfile),
    "small.bk"
  )
  # (gen_tibble() gives a new name instead, see change_duplicated_file_name())
  expect_error(vcf_to_fbm_cpp(vcf_path, quiet = TRUE), "already exists")
  # a backing file can be given with its extension
  other_gt <- vcf_to_fbm_cpp(vcf_path,
    backingfile = file.path(vcf_dir, "other.bk"), quiet = TRUE
  )
  expect_identical(
    basename(attr(other_gt$genotypes, "fbm")$backingfile),
    "other.bk"
  )
})

test_that("the cpp parser finds chromosomes that come back after others", {
  # (the index of each chromosome is looked up in a hash table of the names
  # seen so far)
  chrom <- c("ctg2", "ctg1", "ctg2", "ctg10", "ctg1")
  vcf_path <- tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.3",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb",
      paste(chrom, 1:5, paste0("rs", 1:5), "A\tT\t.\t.\t.\tGT\t0/1\t1/1",
        sep = "\t"
      )
    ),
    vcf_path
  )
  test_gt <- gen_tibble(vcf_path,
    parser = "cpp", backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(
    show_loci(test_gt)$chromosome,
    factor(chrom, levels = c("ctg1", "ctg10", "ctg2"))
  )
})
