test_that("cast_chromosome_to_fac handles diverse cases",{
  #######
  # Basic numeric input
  chr <- c(1,23,2,34,4,3)
  fac_chr <- cast_chromosome_to_factor(chr)
  expect_true(is.factor(fac_chr))
  expect_equal(levels(fac_chr), as.character(c(1,2,3,4,23,34)))
  expect_true(all(fac_chr == chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  expect_true(all(casted_chr == chr))
  #######

  #######
  # Basic character input
  char_chr <- paste0("chr", (chr))
  fac_chr <- cast_chromosome_to_factor(chr)
  expect_true(is.factor(fac_chr))
  expect_equal(levels(fac_chr), as.character(c(1,2,3,4,23,34)))
  # remove chr prefix for comparison
  char_chr <- gsub("chr", "", char_chr)
  expect_true(all(fac_chr == char_chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  casted_chr
  #######

  #######
  # Basic character input with underscore
  chr <- c(1,23,2,34,4,3)
  char_chr <- paste0("chr_", (chr))
  fac_chr <- cast_chromosome_to_factor(chr)
  expect_true(is.factor(fac_chr))
  expect_equal(levels(fac_chr), as.character(c(1,2,3,4,23,34)))
  # remove chr prefix for comparison
  char_chr <- gsub("chr_", "", char_chr)
  expect_true(all(fac_chr == char_chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  casted_chr
  #######

  #######
  # Character input with non-digit chromosomes - transform X and Y to 35 36
  chr <- c(1,23,2,34,4,3,"X","Y")
  fac_chr <- cast_chromosome_to_factor(chr)
  expect_true(is.factor(fac_chr))
  expect_equal(levels(fac_chr), as.character(c(1,2,3,4,23,34,"X","Y")))
  expect_true(all(fac_chr == chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  casted_chr
  #######

  #######
  # Mixed numeric and character with duplicates and input out of order
  chr <- c(3,"Y","Y",1,23,2,2,34,4,"X")
  fac_chr <- cast_chromosome_to_factor(chr)
  expect_true(is.factor(fac_chr))
  expect_equal(levels(fac_chr), as.character(c(1,2,3,4,23,34,"Y","X")))
  expect_true(all(fac_chr == chr))

  # Back to integer
  casted_chr <- cast_chromosome_to_int(fac_chr)
  expect_true(is.integer(casted_chr))
  casted_chr
  #######
})


