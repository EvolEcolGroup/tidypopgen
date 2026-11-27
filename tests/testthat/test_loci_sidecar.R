# delete the backing files of this file when it ends (see helper_cleanup.R)
local_tempdir_cleanup()

anolis_vcf <- system.file(
  "extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
  package = "tidypopgen"
)

# a small VCF: 5 loci on 2 chromosomes ("chr1" gets the code 1, and
# "scaf_2", which is not a number, the code 2)
small_loci <- data.frame(
  chrom = c("chr1", "chr1", "chr1", "scaf_2", "scaf_2"),
  pos = c(10, 20, 30, 5, 15), id = paste0("rs", 1:5), ref = "A", alt = "G"
)
small_vcf <- write_merge_vcf(c("s1", "s2"), small_loci, matrix("0/1", 5, 2))

# write lines to a file with "\n" line ends on every OS (unlike
# writeLines()); with last_eol = "", the last line has no line end
write_raw_lines <- function(path, lines, last_eol = "\n") {
  writeBin(charToRaw(paste0(paste(lines, collapse = "\n"), last_eol)), path)
}

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

test_that("corrupted, truncated or swapped sidecars give clear errors", {
  int_gt <- gen_tibble(small_vcf,
    quiet = TRUE, backingfile = tempfile("small_"),
    names_as_int = TRUE, chromosomes_as_int = TRUE
  )
  path <- loci_sidecar_path(int_gt)
  good <- readLines(path)
  expect_identical(gt_lookup_loci(int_gt, 1:5)$name, paste0("rs", 1:5))
  # the chromosome sidecar copied over the locus sidecar (their headers tell
  # them apart)
  file.copy(chromosomes_sidecar_path(int_gt), path, overwrite = TRUE)
  expect_error(
    gt_lookup_loci(int_gt, 1),
    "not a sidecar file of tidypopgen of the expected type"
  )
  # a copy cut short (its last line has no line end)
  write_raw_lines(path, good, last_eol = "")
  expect_error(gt_lookup_loci(int_gt, 1), "is truncated")
  # a key that is not a number
  bad <- good
  bad[3] <- sub("^2\t", "two\t", bad[3])
  write_raw_lines(path, bad)
  expect_error(
    gt_lookup_loci(int_gt, 1),
    "first field is not a positive integer key"
  )
  # a position that is not a number, which is only read to check positions
  bad <- good
  bad[4] <- sub("\t30$", "\t30x", bad[4])
  write_raw_lines(path, bad)
  expect_error(
    gt_lookup_loci(int_gt, show_loci(int_gt)[3, ]),
    "has a value that is not an integer"
  )
  expect_identical(
    gt_lookup_loci(int_gt, show_loci(int_gt)[3, ],
      check_positions = FALSE
    )$name,
    "rs3"
  )
  # the good sidecar is indexed again, and read as before
  write_raw_lines(path, good)
  expect_identical(gt_lookup_loci(int_gt, 1:5)$name, paste0("rs", 1:5))
})

test_that("sidecars are copied whole for a new backing file, or not at all", {
  int_gt <- gen_tibble(small_vcf,
    quiet = TRUE, backingfile = tempfile("small_"),
    names_as_int = TRUE, chromosomes_as_int = TRUE
  )
  # sidecars left at the new path by older data are replaced
  new_bk <- tempfile("sub_")
  writeLines("stale", paste0(new_bk, ".loci.tsv"))
  writeLines("stale", paste0(new_bk, ".chromosomes.tsv"))
  sub_gt <- gt_update_backingfile(select_loci(int_gt, c(4, 1)),
    backingfile = new_bk, quiet = TRUE
  )
  looked_up <- gt_lookup_loci(sub_gt, 1:2)
  expect_identical(looked_up$name, c("rs4", "rs1"))
  expect_identical(looked_up$chromosome, c("scaf_2", "chr1"))
  # and removed, if the names are not integers
  names_gt <- gen_tibble(small_vcf, quiet = TRUE, backingfile = tempfile())
  new_bk <- tempfile("sub_")
  writeLines("stale", paste0(new_bk, ".loci.tsv"))
  writeLines("stale", paste0(new_bk, ".chromosomes.tsv"))
  gt_update_backingfile(names_gt, backingfile = new_bk, quiet = TRUE)
  expect_false(file.exists(paste0(new_bk, ".loci.tsv")))
  expect_false(file.exists(paste0(new_bk, ".chromosomes.tsv")))
  # sidecars that lack a locus (or a chromosome) that is kept are not copied
  # (a partial sidecar would make some lookups fail much later)
  loci_lines <- readLines(loci_sidecar_path(int_gt))
  write_raw_lines(loci_sidecar_path(int_gt), loci_lines[-2])
  chromosome_lines <- readLines(chromosomes_sidecar_path(int_gt))
  write_raw_lines(chromosomes_sidecar_path(int_gt), chromosome_lines[-2])
  warns <- capture_warnings(
    sub_gt <- gt_update_backingfile(int_gt,
      backingfile = tempfile("sub_"), quiet = TRUE
    )
  )
  expect_match(warns, "does not have all the loci", all = FALSE)
  expect_match(warns, "does not have all the chromosomes", all = FALSE)
  expect_false(file.exists(loci_sidecar_path(sub_gt)))
  expect_false(file.exists(chromosomes_sidecar_path(sub_gt)))
  expect_error(gt_restore_loci(sub_gt), "there is no locus sidecar")
  expect_error(
    gt_lookup_chromosomes(sub_gt, 1L),
    "there is no chromosome sidecar"
  )
})

test_that("integer names and chromosome codes are checked before lookups", {
  names_gt <- gen_tibble(small_vcf, quiet = TRUE, backingfile = tempfile())
  int_gt <- gen_tibble(small_vcf,
    quiet = TRUE, backingfile = tempfile("small_"),
    names_as_int = TRUE, chromosomes_as_int = TRUE
  )
  # names (a factor or strings) are returned as they are, so that the same
  # code works on data read without the integer options
  expect_identical(
    gt_lookup_chromosomes(names_gt, show_loci(names_gt)$chromosome),
    as.character(show_loci(names_gt)$chromosome)
  )
  expect_identical(gt_lookup_chromosomes(int_gt, c("chr1", NA)), c("chr1", NA))
  expect_error(gt_lookup_chromosomes(int_gt, TRUE), "integer codes")
  # codes must be positive whole numbers (doubles are fine)
  expect_identical(
    gt_lookup_chromosomes(int_gt, c(2, 1, 2)),
    c("scaf_2", "chr1", "scaf_2")
  )
  expect_error(gt_lookup_chromosomes(int_gt, 1.5), "should be integers")
  expect_error(gt_lookup_chromosomes(int_gt, 0), "positive integers")
  # a missing name (e.g. from a join that did not find a locus)
  loci <- show_loci(int_gt)[1:2, ]
  loci$name[2] <- NA
  expect_error(gt_lookup_loci(int_gt, loci), "with no NAs")
  # codes replaced by whole numbers stay integer codes, and are restored
  loci <- show_loci(int_gt)
  loci$chromosome <- as.numeric(loci$chromosome)
  show_loci(int_gt) <- loci
  expect_true(is.integer(show_loci(int_gt)$chromosome))
  expect_identical(show_loci(gt_restore_loci(int_gt)), show_loci(names_gt))
  # a code that is not in the chromosome sidecar can not be restored
  loci$chromosome[5] <- 99
  show_loci(int_gt) <- loci
  expect_error(
    gt_restore_loci(int_gt),
    "has no line for 1 of the chromosome codes"
  )
  # .x must be a gen_tibble (or, for lookups, its genotypes)
  expect_error(gt_restore_loci(int_gt$genotypes), "should be a gen_tibble")
  expect_error(gt_lookup_loci(show_loci(int_gt), 1), "should be a gen_tibble")
})

test_that("restoring a subset of loci keeps only their chromosomes", {
  names_gt <- gen_tibble(small_vcf, quiet = TRUE, backingfile = tempfile())
  int_gt <- gen_tibble(small_vcf,
    quiet = TRUE, backingfile = tempfile("small_"),
    names_as_int = TRUE, chromosomes_as_int = TRUE
  )
  # the loci of "scaf_2" only, in reverse order, and with a missing chromosome
  sub_gt <- select_loci(int_gt, c(5, 4))
  restored <- show_loci(gt_restore_loci(sub_gt))
  expected <- show_loci(select_loci(names_gt, c(5, 4)))
  expected$chromosome <- droplevels(expected$chromosome)
  expect_identical(restored, expected)
  expect_identical(levels(restored$chromosome), "scaf_2")
  # missing codes stay missing
  loci <- show_loci(int_gt)
  loci$chromosome[2] <- NA
  show_loci(int_gt) <- loci
  restored <- show_loci(gt_restore_loci(int_gt))$chromosome
  expect_identical(
    as.character(restored),
    c("chr1", NA, "chr1", "scaf_2", "scaf_2")
  )
  expect_identical(levels(restored), c("chr1", "scaf_2"))
})

test_that("chromosomes_as_int is only available for VCF files with cpp", {
  bed_file <- system.file("extdata", "pop_a.bed", package = "tidypopgen")
  expect_error(
    gen_tibble(bed_file, chromosomes_as_int = TRUE, quiet = TRUE),
    "only available for VCF files read with parser = 'cpp'"
  )
  expect_error(
    gen_tibble(small_vcf,
      parser = "vcfR", chromosomes_as_int = TRUE,
      backingfile = tempfile(), quiet = TRUE
    ),
    "only available for VCF files read with parser = 'cpp'"
  )
})

test_that("gt_save lists the sidecars to keep with the files", {
  int_gt <- gen_tibble(small_vcf,
    quiet = TRUE, backingfile = tempfile("small_"),
    names_as_int = TRUE, chromosomes_as_int = TRUE
  )
  msgs <- capture_messages(gt_save(int_gt))
  expect_match(
    msgs, "original names of loci in: .*small_.*\\.loci\\.tsv",
    all = FALSE
  )
  expect_match(
    msgs, "names of chromosomes in: .*small_.*\\.chromosomes\\.tsv",
    all = FALSE
  )
})

test_that("gt_get_file_names and gt_save list the sidecars in use", {
  int_gt <- gen_tibble(small_vcf,
    quiet = TRUE, backingfile = tempfile("small_"),
    names_as_int = TRUE, chromosomes_as_int = TRUE
  )
  files <- gt_get_file_names(int_gt)
  # the .rds and .bk files keep their positions, the sidecars follow (loci
  # first), and all of them exist
  expect_length(files, 4)
  expect_match(files[1], "\\.rds$")
  expect_match(files[2], "\\.bk$")
  expect_identical(
    files[3:4],
    c(loci_sidecar_path(int_gt), chromosomes_sidecar_path(int_gt))
  )
  expect_true(all(file.exists(files)))
  # the genotypes column gives the same files
  expect_identical(gt_get_file_names(int_gt$genotypes), files)
  # gt_save() returns the .gt file, then the same files
  expect_identical(gt_save(int_gt, quiet = TRUE)[-1], files)

  # only the sidecar in use is listed
  name_gt <- gen_tibble(small_vcf,
    quiet = TRUE, backingfile = tempfile("small_"), names_as_int = TRUE
  )
  expect_identical(
    gt_get_file_names(name_gt)[-(1:2)], loci_sidecar_path(name_gt)
  )
  # without integer names or chromosomes, only the two files of the FBM
  plain_gt <- gen_tibble(small_vcf,
    quiet = TRUE, backingfile = tempfile("small_")
  )
  expect_length(gt_get_file_names(plain_gt), 2)
  # a sidecar that is missing is not listed
  unlink(loci_sidecar_path(name_gt))
  expect_length(gt_get_file_names(name_gt), 2)
})

test_that("very long names are kept by the sidecar", {
  # a name of 2 MB, longer than the buffer of the sidecar reader (1 MB)
  long_name <- strrep("x", 2^21)
  loci <- data.frame(
    chrom = "1", pos = 1:3, id = c("rs1", long_name, "rs3"),
    ref = "A", alt = "G"
  )
  vcf <- write_merge_vcf("s1", loci, matrix("0/1", 3, 1))
  int_gt <- gen_tibble(vcf,
    quiet = TRUE, backingfile = tempfile("long_"),
    names_as_int = TRUE
  )
  expect_identical(
    gt_lookup_loci(int_gt, c(3, 2, 1))$name,
    c("rs3", long_name, "rs1")
  )
  expect_identical(
    show_loci(gt_restore_loci(int_gt))$name,
    c("rs1", long_name, "rs3")
  )
})

test_that("missing chromosomes stay missing as integer codes", {
  loci <- data.frame(
    chrom = ".", pos = c(5, 9, 12), id = c("a", "b", "c"),
    ref = "A", alt = "G"
  )
  vcf <- write_merge_vcf(c("s1", "s2"), loci, matrix("0/1", 3, 2))
  names_gt <- gen_tibble(vcf, quiet = TRUE, backingfile = tempfile())
  int_gt <- gen_tibble(vcf,
    quiet = TRUE, backingfile = tempfile(),
    chromosomes_as_int = TRUE
  )
  expect_identical(show_loci(int_gt)$chromosome, rep(NA_integer_, 3))
  # the chromosome sidecar only has its header
  expect_identical(
    readLines(chromosomes_sidecar_path(int_gt)),
    chromosomes_sidecar_header
  )
  expect_identical(show_loci(gt_restore_loci(int_gt)), show_loci(names_gt))
  expect_identical(
    is_loci_table_ordered(int_gt),
    is_loci_table_ordered(names_gt)
  )
  # duplicated positions are found among missing chromosomes too
  loci$pos[2] <- 5
  vcf <- write_merge_vcf(c("s1", "s2"), loci, matrix("0/1", 3, 2))
  for (chromosomes_as_int in c(FALSE, TRUE)) {
    expect_error(
      gen_tibble(vcf,
        quiet = TRUE, backingfile = tempfile(),
        chromosomes_as_int = chromosomes_as_int
      ),
      "duplicated loci"
    )
  }
})

test_that("large chromosome codes are grouped as their names", {
  # chromosome names that are large numbers keep them as their codes (here
  # above 2e8), which are too sparse to be the groups of the C++ helpers
  vcf_lines <- readLines(anolis_vcf)
  body <- which(!startsWith(vcf_lines, "#"))
  vcf_lines[body] <- paste0(
    "chr", 200000000L + as.integer(sub("^locus_([0-9]+)\t.*", "\\1",
      vcf_lines[body]
    )),
    sub("^[^\t]*", "", vcf_lines[body])
  )
  big_vcf <- tempfile(fileext = ".vcf")
  writeLines(vcf_lines, big_vcf)
  names_gt <- gen_tibble(big_vcf, quiet = TRUE, backingfile = tempfile("big_"))
  int_gt <- gen_tibble(big_vcf,
    quiet = TRUE, backingfile = tempfile("big_"),
    chromosomes_as_int = TRUE
  )
  codes <- show_loci(int_gt)$chromosome
  expect_gt(min(codes), 2e8)
  expect_identical(
    codes,
    as.integer(sub("chr", "", show_loci(names_gt)$chromosome))
  )
  # the checks of the order of the loci agree: in order, with the first
  # chromosome split in two runs (it has several loci), and reversed
  n <- count_loci(int_gt)
  ordered <- vapply(
    list(seq_len(n), c(2:n, 1), rev(seq_len(n))),
    function(loci_order) {
      c(
        is_loci_table_ordered(select_loci(int_gt, all_of(loci_order))),
        is_loci_table_ordered(select_loci(names_gt, all_of(loci_order)))
      )
    },
    logical(2)
  )
  expect_identical(ordered[1, ], ordered[2, ])
  expect_identical(ordered[1, ], c(TRUE, FALSE, FALSE))
  # runs of homozygosity are the same, with the names of the chromosomes
  # (from the chromosome sidecar)
  roh_args <- list(
    window_snp = 2, window_het = 0, window_missing = 2, min_snp = 2,
    min_kb = 0, max_kb_per_snp = 1000, n_cores = 1
  )
  names_roh <- do.call(windows_indiv_roh, c(list(names_gt), roh_args))
  int_roh <- do.call(windows_indiv_roh, c(list(int_gt), roh_args))
  expect_gt(nrow(int_roh), 0)
  expect_identical(int_roh, names_roh)
  # without the sidecar, the runs have the codes (sorted as the names here,
  # as all the names have the same number of digits)
  unlink(chromosomes_sidecar_path(int_gt))
  int_roh <- do.call(windows_indiv_roh, c(list(int_gt), roh_args))
  expect_true(all(int_roh$chromosome > 2e8))
  expect_identical(
    int_roh$chromosome,
    as.integer(sub("chr", "", names_roh$chromosome))
  )
  # (the attribute chromosome_length has the codes too)
  expect_equal(
    dplyr::select(int_roh, -"chromosome"),
    dplyr::select(names_roh, -"chromosome"),
    ignore_attr = TRUE
  )
  # a duplicated locus (the same chromosome and position) is found
  dup <- strsplit(vcf_lines[body[1]], "\t")[[1]]
  dup[3] <- "dup_id"
  dup_lines <- append(vcf_lines, paste(dup, collapse = "\t"), after = body[1])
  writeLines(dup_lines, big_vcf)
  for (chromosomes_as_int in c(FALSE, TRUE)) {
    expect_error(
      gen_tibble(big_vcf,
        quiet = TRUE, backingfile = tempfile(),
        chromosomes_as_int = chromosomes_as_int
      ),
      "duplicated loci"
    )
  }
})

test_that("sidecar copies and reads check the order and number of lines", {
  int_gt <- gen_tibble(small_vcf,
    quiet = TRUE, backingfile = tempfile("small_"),
    names_as_int = TRUE
  )
  path <- loci_sidecar_path(int_gt)
  # whole columns, counting the lines (or checking their number, if given)
  read_cols <- function(...) {
    sidecar_read(path, loci_sidecar_header,
      fields = c(1L, 2L), types = c("integer", "character"), ...
    )
  }
  cols <- read_cols()
  expect_identical(cols, list(1:5, paste0("rs", 1:5)))
  expect_identical(read_cols(n = 5), cols)
  expect_error(read_cols(n = 4), "more lines than expected")
  expect_error(read_cols(n = 6), "fewer lines than expected")
  # a file whose keys are not in order is never copied (and the partial copy
  # is removed)
  unsorted <- tempfile()
  write_raw_lines(unsorted, readLines(path)[c(1, 3, 2, 4:6)])
  out <- tempfile()
  expect_error(
    sidecar_subset(unsorted, out, loci_sidecar_header, 1:2),
    "not strictly increasing"
  )
  expect_false(file.exists(out))
  expect_error(
    sidecar_gather(unsorted, loci_sidecar_header,
      index_keys = integer(0), index_offsets = numeric(0),
      keys = integer(0), path_out = out, first_new_key = 1, append = FALSE
    ),
    "not strictly increasing"
  )
  expect_false(file.exists(out))
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
  # names that are only a prefix are other names, and NA is missing (as ".")
  expect_identical(
    chromosome_int_codes_cpp(c("chr", "chr3", NA, "Chromosome_")),
    c(4L, 3L, NA, 5L)
  )
  # the largest integer can be a code, but then no other name can be coded
  expect_identical(
    chromosome_int_codes_cpp(c("2147483646", "chrX")),
    c(2147483646L, 2147483647L)
  )
  expect_error(
    chromosome_int_codes_cpp(c("2147483647", "chrX")),
    "too many chromosomes"
  )
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
  # (replacement codes must be whole numbers, rather than being truncated)
  new_loci <- show_loci(int_gt)
  new_loci$chromosome <- new_loci$chromosome + 0.5
  expect_error(show_loci(int_gt) <- new_loci, "whole numbers")
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
  # the same runs, in the same order, with the names of the chromosomes
  # (also in the attribute with the lengths of the chromosomes)
  expect_identical(int_roh, names_roh)
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
  # (a key repeated within a chunk)
  out_repeated <- tempfile()
  sidecar_gather(path, loci_sidecar_header,
    index_keys = index$keys, index_offsets = index$offsets,
    keys = c(3L, 3L, 1L), path_out = out_repeated, first_new_key = 1,
    append = FALSE
  )
  expect_identical(readLines(out_repeated)[-1], c(
    "1\trs3\t1\t3", "2\trs3\t1\t3", "3\trs1\t1\t1"
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
