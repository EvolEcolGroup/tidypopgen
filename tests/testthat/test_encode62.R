test_that("encode62 and decode62 work correctly", {
  test_chrom <- c(1,10,260)
  test_pos <- c(1,1000,1000000)
  
  encoded_values <- encode62(test_chrom, test_pos)
  decoded_values <- decode62(encoded_values, max_chr = 260)
  expect_equal(decoded_values[,1], test_chrom)
  expect_equal(decoded_values[,2], test_pos)
})
