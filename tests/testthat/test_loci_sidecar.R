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
    "should have a 'name' column"
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

test_that("sidecars left by older files with the same name are removed", {
  backingfile <- tempfile("anolis_")
  stale <- paste0(backingfile, ".loci.tsv")
  writeLines("stale", stale)
  names_gt <- gen_tibble(anolis_vcf, quiet = TRUE, backingfile = backingfile)
  expect_false(file.exists(stale))
})
