test_that("cast_chromosome_to_factor handles diverse cases", {
  #######
  # Basic numeric input
  chr <- c(1, 23, 2, 34, 4, 3)
  fac_chr <- cast_chromosome_to_factor(chr)
  expect_true(is.factor(fac_chr))
  expect_equal(levels(fac_chr), as.character(c(1, 2, 3, 4, 23, 34)))
  expect_equal(as.character(fac_chr), as.character(chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  expect_equal(casted_chr, chr)
  #######

  #######
  # Basic character input
  char_chr <- paste0("chr", (chr))
  fac_chr <- cast_chromosome_to_factor(char_chr)
  expect_true(is.factor(fac_chr))
  expect_equal(
    levels(fac_chr),
    c("chr1", "chr2", "chr3", "chr4", "chr23", "chr34")
  )
  expect_true(all(fac_chr == char_chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  expect_equal(casted_chr, chr)
  #######

  #######
  # Basic character input with underscore
  chr <- c(1, 23, 2, 34, 4, 3)
  char_chr <- paste0("chr_", (chr))
  fac_chr <- cast_chromosome_to_factor(char_chr)
  expect_true(is.factor(fac_chr))
  expect_equal(
    levels(fac_chr),
    c("chr_1", "chr_2", "chr_3", "chr_4", "chr_23", "chr_34")
  )
  expect_true(all(fac_chr == char_chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  expect_equal(casted_chr, chr)
  #######

  #######
  # Character input with non-digit chromosomes - transform X and Y to 35 36
  chr <- c(1, 23, 2, 34, 4, 3, "X", "Y")
  fac_chr <- cast_chromosome_to_factor(chr)
  expect_true(is.factor(fac_chr))
  expect_equal(levels(fac_chr), as.character(c(1, 2, 3, 4, 23, 34, "X", "Y")))
  expect_true(all(fac_chr == chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  expect_equal(casted_chr, c(1, 23, 2, 34, 4, 3, 35, 36))
  #######

  #######
  # Mixed numeric and character with duplicates and input out of order
  chr <- c(3, "Y", "Y", 1, 23, 2, 2, 34, 4, "X")
  fac_chr <- cast_chromosome_to_factor(chr)
  expect_true(is.factor(fac_chr))
  expect_equal(levels(fac_chr), as.character(c(1, 2, 3, 4, 23, 34, "X", "Y")))
  expect_true(all(fac_chr == chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  expect_equal(casted_chr, c(3, 36, 36, 1, 23, 2, 2, 34, 4, 35))
  #######

  #######
  # Mixed numeric and character with duplicates and input out of order
  chr <- c(
    "chr3", "chrY", "chrY",
    "chr1", "chr23", "chr2",
    "chr2", "chr34", "chr4", "chrX"
  )
  fac_chr <- cast_chromosome_to_factor(chr)
  expect_true(is.factor(fac_chr))
  expect_equal(
    levels(fac_chr),
    c("chr1", "chr2", "chr3", "chr4", "chr23", "chr34", "chrX", "chrY")
  )
  expect_true(all(fac_chr == chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  expect_equal(casted_chr, c(3, 36, 36, 1, 23, 2, 2, 34, 4, 35))
  #######

  #######
  # Inputs with NA values
  chr <- c(3, NA, "Y", 1, 23, NA, 2, 34, 4, "X")
  expect_error(
    cast_chromosome_to_factor(chr),
    "NA values are not allowed in chromosome names"
  )
  #######

  #######
  # Only non-digits
  chr <- c("X", "Y", "MT")
  fac_chr <- cast_chromosome_to_factor(chr)
  expect_true(is.factor(fac_chr))
  expect_equal(levels(fac_chr), as.character(c("MT", "X", "Y")))
  expect_true(all(fac_chr == chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  expect_equal(casted_chr, c(2, 3, 1))
  #######
})
