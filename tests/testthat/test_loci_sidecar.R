anolis_vcf <- system.file(
  "extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
  package = "tidypopgen"
)

test_that("the locus sidecar restores the names of names_as_int = FALSE", {
  names_gt <- gen_tibble(anolis_vcf,
    quiet = TRUE, backingfile = tempfile("anolis_")
  )
  # names are not integers: no sidecar
  expect_false(file.exists(loci_sidecar_path(names_gt)))
  for (n_cores in c(1, 2)) {
    int_gt <- gen_tibble(anolis_vcf,
      quiet = TRUE, backingfile = tempfile("anolis_"),
      names_as_int = TRUE, n_cores = n_cores
    )
    expect_true(file.exists(loci_sidecar_path(int_gt)))
    expect_identical(
      show_loci(gt_restore_loci(int_gt)),
      show_loci(names_gt)
    )
  }
  # the sidecar can be read as a plain table
  sidecar <- utils::read.delim(loci_sidecar_path(int_gt), quote = "")
  expect_identical(names(sidecar), c(
    "int_name", "name", "chromosome", "position"
  ))
  expect_identical(sidecar$int_name, show_loci(int_gt)$name)
  expect_identical(sidecar$name, show_loci(names_gt)$name)
})

test_that("missing IDs and chromosomes are restored as names_as_int = FALSE", {
  # missing IDs (filled with the chromosome and position, NA for a missing
  # chromosome), a missing chromosome, and an indel (not read)
  loci <- data.frame(
    chrom = c("chr1", "chr1", "chr1", ".", "chr2"),
    pos = c(5, 10, 12, 7, 3),
    id = c("rs1", ".", "rs3", ".", "."),
    ref = c("A", "C", "AT", "G", "T"),
    alt = c("G", "T", "A", "A", "C")
  )
  genotypes <- matrix(c("0/1", "1/1"), nrow = 5, ncol = 2, byrow = TRUE)
  vcf <- write_merge_vcf(c("s1", "s2"), loci, genotypes)
  names_gt <- gen_tibble(vcf, quiet = TRUE, backingfile = tempfile("vcf_"))
  int_gt <- gen_tibble(vcf,
    quiet = TRUE, backingfile = tempfile("vcf_"),
    names_as_int = TRUE
  )
  expect_identical(
    show_loci(names_gt)$name,
    c("rs1", "chr1_10", "NA_7", "chr2_3")
  )
  expect_identical(
    show_loci(gt_restore_loci(int_gt)),
    show_loci(names_gt)
  )
  # a missing chromosome is "." in the sidecar
  expect_identical(
    readLines(loci_sidecar_path(int_gt))[4],
    "3\tNA_7\t.\t7"
  )
})

test_that("gt_lookup_loci looks up rows and subsets of the loci table", {
  # (a tiny block size, so that lookups cross the edges of blocks)
  withr::local_options(tidypopgen.sidecar_block = 7L)
  names_gt <- gen_tibble(anolis_vcf,
    quiet = TRUE, backingfile = tempfile("anolis_")
  )
  int_gt <- gen_tibble(anolis_vcf,
    quiet = TRUE, backingfile = tempfile("anolis_"),
    names_as_int = TRUE
  )
  all_names <- show_loci(names_gt)$name
  # row numbers, unsorted and repeated
  rows <- c(50, 3, 3, count_loci(int_gt), 1, 29)
  expect_identical(gt_lookup_loci(int_gt, rows)$name, all_names[rows])
  # the genotypes column works too
  expect_identical(
    gt_lookup_loci(int_gt$genotypes, rows)$name,
    all_names[rows]
  )
  # a subset of the loci table, in another order, with an extra column, is
  # the same table with the names translated
  hits <- show_loci(int_gt)
  hits$stat <- seq_len(nrow(hits))
  hits <- hits[rev(rows), ]
  res <- gt_lookup_loci(int_gt, hits)
  expect_identical(res$name, all_names[rev(rows)])
  expect_identical(res[, -2], hits[, -2])
  # the row form equals the table form
  expect_identical(
    gt_lookup_loci(int_gt, rows),
    gt_lookup_loci(int_gt, show_loci(int_gt)[rows, ])
  )
  # tables whose names are not integers are returned unchanged
  expect_identical(
    gt_lookup_loci(names_gt, show_loci(names_gt)[rows, ]),
    show_loci(names_gt)[rows, ]
  )
  expect_identical(gt_restore_loci(names_gt), names_gt)
  # errors
  expect_error(gt_lookup_loci(int_gt, 0), "row numbers")
  expect_error(gt_lookup_loci(int_gt, "a"), "subset of the loci table")
  expect_error(
    gt_lookup_loci(int_gt, data.frame(x = 1)),
    "should have a .name. or .chromosome. column"
  )
})

test_that("the locus sidecar follows subsets, reordering and rbind", {
  names_gt <- gen_tibble(anolis_vcf,
    quiet = TRUE, backingfile = tempfile("anolis_")
  )
  int_gt <- gen_tibble(anolis_vcf,
    quiet = TRUE, backingfile = tempfile("anolis_"),
    names_as_int = TRUE
  )
  # select and reorder loci, then write a new backing file
  keep <- sample(count_loci(int_gt), 100)
  sub_gt <- int_gt %>%
    select_loci(all_of(keep)) %>%
    gt_update_backingfile(backingfile = tempfile("anolis_sub_"), quiet = TRUE)
  # the new sidecar only has the loci that are kept, in order of their keys
  expect_identical(
    length(readLines(loci_sidecar_path(sub_gt))),
    101L
  )
  expect_identical(
    show_loci(gt_restore_loci(sub_gt))$name,
    show_loci(names_gt)$name[keep]
  )
  # gt_order_loci() (which writes a new backing file)
  ordered_gt <- gt_order_loci(sub_gt, use_current_table = FALSE, quiet = TRUE)
  expect_identical(
    gt_lookup_loci(ordered_gt, show_loci(ordered_gt))$name,
    show_loci(gt_order_loci(names_gt %>% select_loci(all_of(keep)),
      use_current_table = FALSE, quiet = TRUE
    ))$name
  )

  # rbind keeps the loci of ref (matched by position with integer names)
  pop_a_int <- gen_tibble(
    system.file("extdata/pop_a.vcf", package = "tidypopgen"),
    quiet = TRUE, backingfile = tempfile("pop_a_"), names_as_int = TRUE
  )
  pop_b_int <- gen_tibble(
    system.file("extdata/pop_b.vcf", package = "tidypopgen"),
    quiet = TRUE, backingfile = tempfile("pop_b_"), names_as_int = TRUE
  )
  pop_a <- gen_tibble(
    system.file("extdata/pop_a.vcf", package = "tidypopgen"),
    quiet = TRUE, backingfile = tempfile("pop_a_")
  )
  merged_gt <- rbind(pop_a_int, pop_b_int,
    flip_strand = TRUE, quiet = TRUE,
    backingfile = tempfile("merged_")
  )
  expect_true(file.exists(loci_sidecar_path(merged_gt)))
  merged_loci <- show_loci(gt_restore_loci(merged_gt))
  pop_a_loci <- show_loci(pop_a)
  expect_identical(
    merged_loci$name,
    pop_a_loci$name[match(merged_loci$position, pop_a_loci$position)]
  )
})

test_that("lookups check that the sidecar belongs to the data", {
  int_gt <- gen_tibble(anolis_vcf,
    quiet = TRUE, backingfile = tempfile("anolis_"),
    names_as_int = TRUE
  )
  # a different position is an error, unless check_positions = FALSE
  loci <- show_loci(int_gt)[1:5, ]
  loci$position[3] <- loci$position[3] + 1L
  expect_error(gt_lookup_loci(int_gt, loci), "does not match the loci")
  expect_identical(
    gt_lookup_loci(int_gt, loci, check_positions = FALSE)$name,
    gt_lookup_loci(int_gt, 1:5)$name
  )
  # so is a locus that is not in the sidecar
  loci <- show_loci(int_gt)[1:2, ]
  loci$name[2] <- 10L^6
  expect_error(gt_lookup_loci(int_gt, loci), "has no line for 1 of the loci")
  # a sidecar edited by hand (with unsorted keys) is an error
  path <- loci_sidecar_path(int_gt)
  lines <- readLines(path)
  writeLines(lines[c(1, 3, 2, 4:length(lines))], path)
  expect_error(gt_lookup_loci(int_gt, 1), "not strictly increasing")
  # as is a missing sidecar
  unlink(path)
  expect_error(gt_restore_loci(int_gt), "there is no locus sidecar")
  # a new backing file then has no sidecar (and lookups give the error)
  sub_gt <- gt_update_backingfile(int_gt,
    backingfile = tempfile("anolis_sub_"), quiet = TRUE
  )
  expect_false(file.exists(loci_sidecar_path(sub_gt)))
})

test_that("sidecars with Windows line endings can be read", {
  int_gt <- gen_tibble(anolis_vcf,
    quiet = TRUE, backingfile = tempfile("anolis_"),
    names_as_int = TRUE
  )
  expected <- gt_lookup_loci(int_gt, 1:5)
  path <- loci_sidecar_path(int_gt)
  lines <- readLines(path)
  con <- file(path, open = "wb")
  writeLines(lines, con, sep = "\r\n")
  close(con)
  expect_identical(gt_lookup_loci(int_gt, 1:5), expected)
})

test_that("sidecars left by older files with the same name are removed", {
  backingfile <- tempfile("anolis_")
  stale <- paste0(backingfile, ".loci.tsv")
  writeLines("stale", stale)
  names_gt <- gen_tibble(anolis_vcf, quiet = TRUE, backingfile = backingfile)
  expect_false(file.exists(stale))
})

test_that("chromosome names that are numbers keep their value as codes", {
  chromosomes <- c(
    "chr2", "1", "chr10", "chrX", "chr1", ".", "scaffold_5", "0",
    "Chromosome_7", "99999999999"
  )
  # numbers keep their value (with a chr or chromosome prefix); the second
  # name with the same number, 0, numbers too large for an integer, and other
  # names are coded after the largest number, in order of appearance; "." is
  # missing
  expect_identical(
    chromosome_int_codes_cpp(chromosomes),
    c(2L, 1L, 10L, 11L, 12L, NA, 13L, 14L, 7L, 15L)
  )
  # the same codes when reading a VCF
  loci <- data.frame(
    chrom = chromosomes,
    pos = seq_along(chromosomes),
    id = paste0("rs", seq_along(chromosomes)),
    ref = "A", alt = "G"
  )
  genotypes <- matrix("0/1", nrow = length(chromosomes), ncol = 2)
  vcf <- write_merge_vcf(c("s1", "s2"), loci, genotypes)
  names_gt <- gen_tibble(vcf, quiet = TRUE, backingfile = tempfile("vcf_"))
  int_gt <- gen_tibble(vcf,
    quiet = TRUE, backingfile = tempfile("vcf_"),
    chromosomes_as_int = TRUE
  )
  expect_identical(
    show_loci(int_gt)$chromosome,
    chromosome_int_codes_cpp(chromosomes)
  )
  # the chromosome sidecar is sorted by code, without the missing chromosome
  sidecar <- utils::read.delim(chromosomes_sidecar_path(int_gt), quote = "")
  expect_identical(sidecar$int_chromosome, sort(na.omit(
    chromosome_int_codes_cpp(chromosomes)
  ))[seq_len(nrow(sidecar))])
  expect_identical(
    show_loci(gt_restore_loci(int_gt)),
    show_loci(names_gt)
  )
  expect_identical(
    gt_lookup_chromosomes(int_gt, c(14L, NA, 2L)),
    c("0", NA, "chr2")
  )
  expect_error(gt_lookup_chromosomes(int_gt, 3L), "has no line")
})

test_that("integer chromosomes are restored, looked up and kept", {
  names_gt <- gen_tibble(anolis_vcf,
    quiet = TRUE, backingfile = tempfile("anolis_")
  )
  for (names_as_int in c(FALSE, TRUE)) {
    int_gt <- gen_tibble(anolis_vcf,
      quiet = TRUE, backingfile = tempfile("anolis_"),
      names_as_int = names_as_int, chromosomes_as_int = TRUE
    )
    expect_true(is.integer(show_loci(int_gt)$chromosome))
    expect_identical(
      show_loci(gt_restore_loci(int_gt)),
      show_loci(names_gt)
    )
    # only the chromosomes, or only the names
    expect_true(is.integer(show_loci(
      gt_restore_loci(int_gt, chromosomes = FALSE)
    )$chromosome))
    expect_identical(
      show_loci(gt_restore_loci(int_gt, names = FALSE))$chromosome,
      show_loci(names_gt)$chromosome
    )
    rows <- c(30, 2, 2, 700)
    expect_identical(
      gt_lookup_loci(int_gt, rows)$chromosome,
      as.character(show_loci(names_gt)$chromosome[rows])
    )
    expect_identical(
      gt_lookup_loci(int_gt, rows)$name,
      show_loci(names_gt)$name[rows]
    )
  }
  # the integer chromosomes are kept by the checks and functions of the loci
  # table
  expect_true(is_loci_table_ordered(int_gt))
  keep <- sample(count_loci(int_gt), 200)
  sub_gt <- int_gt %>%
    select_loci(all_of(keep)) %>%
    gt_order_loci(use_current_table = FALSE, quiet = TRUE)
  expect_true(is.integer(show_loci(sub_gt)$chromosome))
  # the new chromosome sidecar has only the chromosomes that are kept
  n_chromosomes <- length(unique(show_loci(sub_gt)$chromosome))
  expect_identical(
    length(readLines(chromosomes_sidecar_path(sub_gt))),
    n_chromosomes + 1L
  )
  restored <- show_loci(gt_restore_loci(sub_gt))
  names_loci <- show_loci(names_gt)
  names_loci <- names_loci[match(restored$name, names_loci$name), ]
  expect_identical(
    as.character(restored$chromosome),
    as.character(names_loci$chromosome)
  )
})

test_that("rbind matches integer chromosomes by name", {
  pop_a <- system.file("extdata/pop_a.vcf", package = "tidypopgen")
  pop_b <- system.file("extdata/pop_b.vcf", package = "tidypopgen")
  ref_names <- gen_tibble(pop_a, quiet = TRUE, backingfile = tempfile("a_"))
  target_names <- gen_tibble(pop_b,
    quiet = TRUE, backingfile = tempfile("b_")
  )
  merged_names <- rbind(ref_names, target_names,
    flip_strand = TRUE, quiet = TRUE, backingfile = tempfile("merged_")
  )
  ref_int <- gen_tibble(pop_a,
    quiet = TRUE, backingfile = tempfile("a_"),
    chromosomes_as_int = TRUE
  )
  target_int <- gen_tibble(pop_b,
    quiet = TRUE, backingfile = tempfile("b_"),
    chromosomes_as_int = TRUE
  )
  # integer chromosomes in both, or in one of the two datasets
  for (target in list(target_int, target_names)) {
    merged_int <- rbind(ref_int, target,
      flip_strand = TRUE, quiet = TRUE, backingfile = tempfile("merged_")
    )
    expect_true(is.integer(show_loci(merged_int)$chromosome))
    expect_identical(
      show_loci(gt_restore_loci(merged_int))$name,
      show_loci(merged_names)$name
    )
    expect_identical(
      as.character(show_loci(gt_restore_loci(merged_int))$chromosome),
      as.character(show_loci(merged_names)$chromosome)
    )
    expect_identical(
      show_genotypes(merged_int),
      show_genotypes(merged_names)
    )
  }
})

test_that("window statistics and runs keep integer chromosomes", {
  names_gt <- gen_tibble(anolis_vcf,
    quiet = TRUE, backingfile = tempfile("anolis_")
  )
  int_gt <- gen_tibble(anolis_vcf,
    quiet = TRUE, backingfile = tempfile("anolis_"),
    chromosomes_as_int = TRUE
  )
  pops <- rep(c("a", "b"), length.out = nrow(names_gt))
  names_fst <- names_gt %>%
    dplyr::mutate(population = pops) %>%
    dplyr::group_by(population) %>%
    windows_pairwise_pop_fst(
      window_size = 3, step_size = 2, size_unit = "snp", min_loci = 2
    )
  int_fst <- int_gt %>%
    dplyr::mutate(population = pops) %>%
    dplyr::group_by(population) %>%
    windows_pairwise_pop_fst(
      window_size = 3, step_size = 2, size_unit = "snp", min_loci = 2
    )
  expect_true(is.integer(int_fst$chromosome))
  expect_equal(
    dplyr::select(int_fst, -"chromosome"),
    dplyr::select(names_fst, -"chromosome")
  )
  expect_identical(
    gt_lookup_chromosomes(int_gt, int_fst$chromosome),
    as.character(names_fst$chromosome)
  )
  # runs of homozygosity (short runs, as RAD loci are short)
  roh_args <- list(
    window_snp = 2, window_het = 0, window_missing = 2, min_snp = 2,
    min_kb = 0, max_kb_per_snp = 1000, n_cores = 1
  )
  names_roh <- do.call(windows_indiv_roh, c(list(names_gt), roh_args))
  int_roh <- do.call(windows_indiv_roh, c(list(int_gt), roh_args))
  expect_gt(nrow(int_roh), 0)
  expect_true(is.integer(int_roh$chromosome))
  # the same runs (ordered by the codes of the chromosomes, rather than by
  # their names)
  int_roh$chromosome <- gt_lookup_chromosomes(int_gt, int_roh$chromosome)
  names_roh$chromosome <- as.character(names_roh$chromosome)
  # (the lengths of the chromosomes, an attribute, have their codes)
  expect_true(is.integer(attr(int_roh, "chromosome_length")$chromosome))
  sort_runs <- function(x) {
    attr(x, "chromosome_length") <- NULL
    as.data.frame(dplyr::arrange(x, .data$id, .data$chromosome, .data$start))
  }
  expect_identical(sort_runs(int_roh), sort_runs(names_roh))
})

# three files with different samples, with missing IDs, several chromosomes
# (in different orders), and sites first found in each file
sidecar_merge_files <- c(
  write_merge_vcf("a",
    data.frame(
      chrom = c("chr2", "chr2", "chr1", "scaf_9"), pos = c(10, 20, 5, 1),
      id = c(".", "a2", "a3", "a4"), ref = "A", alt = "C"
    ),
    genotypes = matrix("0/1", 4, 1)
  ),
  write_merge_vcf("b",
    data.frame(
      chrom = c("chr1", "chr1", "chr2", "chrX", "chr2"),
      pos = c(5, 7, 15, 3, 20),
      id = c("b1", "b2", ".", "b4", "b5"), ref = "A", alt = "C"
    ),
    genotypes = matrix("1/1", 5, 1)
  ),
  write_merge_vcf("c",
    data.frame(
      chrom = c("chr2", "chrX", "chr1", "chr1"), pos = c(20, 3, 1, 5),
      id = c("c1", "c2", "c3", "c4"), ref = "A", alt = "C"
    ),
    genotypes = matrix("0/0", 4, 1)
  )
)

test_that("gt_merge_vcfs writes sidecars that restore the names", {
  withr::local_options(tidypopgen.sidecar_block = 2L)
  files <- sidecar_merge_files
  for (variants in c("union", "intersection")) {
    names_gt <- gt_merge_vcfs(files,
      variants = variants, names_as_int = FALSE,
      backingfile = tempfile(), quiet = TRUE
    )
    for (chromosomes_as_int in c(FALSE, TRUE)) {
      int_gt <- gt_merge_vcfs(files,
        variants = variants, chromosomes_as_int = chromosomes_as_int,
        backingfile = tempfile(), quiet = TRUE
      )
      expect_true(is.integer(show_loci(int_gt)$name))
      expect_identical(
        is.integer(show_loci(int_gt)$chromosome),
        chromosomes_as_int
      )
      # no temporary files are left
      expect_identical(
        list.files(dirname(loci_sidecar_path(int_gt)),
          pattern = basename(paste0(loci_sidecar_path(int_gt), "."))
        ),
        character(0)
      )
      expect_identical(
        show_loci(gt_restore_loci(int_gt)),
        show_loci(names_gt)
      )
      expect_identical(show_genotypes(int_gt), show_genotypes(names_gt))
    }
  }
  # each locus takes the name it has in the first file where it is found
  expect_identical(
    show_loci(names_gt)$name,
    c("a2", "a3")
  )
})

test_that("gt_concat_vcfs writes sidecars that restore the names", {
  withr::local_options(tidypopgen.sidecar_block = 2L)
  files <- c(
    write_merge_vcf(c("a", "b"),
      data.frame(
        chrom = c("1", "1", "scaf_2"), pos = c(100, 200, 5),
        id = c("s1", ".", "s3"), ref = "A", alt = "G"
      ),
      genotypes = matrix("0/1", 3, 2)
    ),
    write_merge_vcf(c("b", "a"),
      data.frame(
        chrom = c("2", ".", "chrX"), pos = c(50, 60, 70),
        id = c("s4", ".", "s6"), ref = "T", alt = "C"
      ),
      genotypes = matrix("1/1", 3, 2)
    )
  )
  names_gt <- gt_concat_vcfs(files,
    names_as_int = FALSE, backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(
    show_loci(names_gt)$name,
    c("s1", "1_200", "s3", "s4", "NA_60", "s6")
  )
  for (chromosomes_as_int in c(FALSE, TRUE)) {
    int_gt <- gt_concat_vcfs(files,
      chromosomes_as_int = chromosomes_as_int,
      backingfile = tempfile(), quiet = TRUE
    )
    expect_identical(
      show_loci(gt_restore_loci(int_gt)),
      show_loci(names_gt)
    )
    expect_false(file.exists(paste0(loci_sidecar_path(int_gt), ".scan")))
  }
  expect_identical(
    show_loci(int_gt)$chromosome,
    c(1L, 1L, 3L, 2L, NA, 4L)
  )
})

test_that("sidecar_gather writes lines in any order, with new keys", {
  loci <- data.frame(
    chrom = "1", pos = 1:5, id = paste0("rs", 1:5), ref = "A", alt = "G"
  )
  vcf <- write_merge_vcf("a", loci, matrix("0/1", 5, 1))
  int_gt <- gen_tibble(vcf,
    names_as_int = TRUE, backingfile = tempfile(), quiet = TRUE
  )
  path <- loci_sidecar_path(int_gt)
  index <- sidecar_index(path, loci_sidecar_header, 2L)
  out <- tempfile()
  # (in chunks of 2 keys, with a repeated key)
  sidecar_gather(path, loci_sidecar_header,
    index_keys = index$keys, index_offsets = index$offsets,
    keys = c(4L, 1L, 4L, 5L, 2L), path_out = out, first_new_key = 1,
    append = FALSE, chunk_size = 2L
  )
  expect_identical(readLines(out), c(
    loci_sidecar_header,
    "1\trs4\t1\t4", "2\trs1\t1\t1", "3\trs4\t1\t4", "4\trs5\t1\t5",
    "5\trs2\t1\t2"
  ))
  # all the lines, appended with keys after the last one
  sidecar_gather(path, loci_sidecar_header,
    index_keys = integer(0), index_offsets = numeric(0),
    keys = integer(0), path_out = out, first_new_key = 6, append = TRUE
  )
  expect_identical(readLines(out)[7:11], paste0(
    6:10, "\trs", 1:5, "\t1\t", 1:5
  ))
  expect_error(
    sidecar_gather(path, loci_sidecar_header,
      index_keys = index$keys, index_offsets = index$offsets,
      keys = 9L, path_out = tempfile(), first_new_key = 1, append = FALSE
    ),
    "has no line for the key 9"
  )
})
