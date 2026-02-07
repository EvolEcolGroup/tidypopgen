test_that("matrix_int_names creation works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c(100L, 200L, 300L, 400L)
  )

  expect_s3_class(im, "matrix_int_names")
  expect_s3_class(im, "matrix")
  expect_equal(dim(im), c(3, 4))
  expect_equal(row_names(im), c(10L, 20L, 30L))
  expect_equal(col_names(im), c(100L, 200L, 300L, 400L))
})

test_that("matrix_int_names requires integer or character names", {
  m <- matrix(1:12, nrow = 3, ncol = 4)

  expect_error(
    matrix_int_names(m, row_names = c(1, 2, 3)),
    "row_names must be integer or character vector"
  )

  expect_error(
    matrix_int_names(m, col_names = list(1, 2, 3, 4)),
    "col_names must be integer or character vector"
  )
})

test_that("matrix_int_names validates name lengths", {
  m <- matrix(1:12, nrow = 3, ncol = 4)

  expect_error(
    matrix_int_names(m, row_names = c(1L, 2L)),
    "Length of row_names must match number of rows"
  )

  expect_error(
    matrix_int_names(m, col_names = c(1L, 2L, 3L)),
    "Length of col_names must match number of columns"
  )
})

test_that("matrix_int_names works without names", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m)

  expect_s3_class(im, "matrix_int_names")
  expect_null(row_names(im))
  expect_null(col_names(im))
})

test_that("subsetting by position works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c(100L, 200L, 300L, 400L)
  )

  # Single element
  expect_equal(im[1, 1], 1)
  expect_equal(im[2, 3], 8)

  # Row subsetting
  result <- im[1:2, ]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(2, 4))
  expect_equal(row_names(result), c(10L, 20L))
  expect_equal(col_names(result), c(100L, 200L, 300L, 400L))

  # Column subsetting
  result <- im[, 2:3]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(3, 2))
  expect_equal(row_names(result), c(10L, 20L, 30L))
  expect_equal(col_names(result), c(200L, 300L))
})

test_that("subsetting by integer names works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c(100L, 200L, 300L, 400L)
  )

  # Subset by row names
  result <- im[, , i_names = c(10L, 30L)]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(2, 4))
  expect_equal(row_names(result), c(10L, 30L))
  expect_equal(as.vector(result[1, ]), c(1, 4, 7, 10))
  expect_equal(as.vector(result[2, ]), c(3, 6, 9, 12))

  # Subset by column names
  result <- im[, , j_names = c(200L, 400L)]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(3, 2))
  expect_equal(col_names(result), c(200L, 400L))
  expect_equal(as.vector(result[, 1]), c(4, 5, 6))
  expect_equal(as.vector(result[, 2]), c(10, 11, 12))

  # Subset by both
  result <- im[, , i_names = c(20L, 30L), j_names = c(100L, 300L)]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(2, 2))
  expect_equal(row_names(result), c(20L, 30L))
  expect_equal(col_names(result), c(100L, 300L))
})

test_that("subsetting with drop works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c(100L, 200L, 300L, 400L)
  )

  # Single row drops to vector by default
  result <- im[1, ]
  expect_true(is.vector(result))
  expect_equal(length(result), 4)

  # Single column drops to vector by default
  result <- im[, 1]
  expect_true(is.vector(result))
  expect_equal(length(result), 3)

  # drop = FALSE preserves matrix
  result <- im[1, , drop = FALSE]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(1, 4))
})

test_that("assignment to matrix_int_names works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c(100L, 200L, 300L, 400L)
  )

  # Assign single element
  im[1, 1] <- 99
  expect_equal(im[1, 1], 99)

  # Assign row
  im[2, ] <- c(88, 77, 66, 55)
  expect_equal(as.vector(im[2, ]), c(88, 77, 66, 55))

  # Assign using integer names
  im[, , i_names = 30L, j_names = 100L] <- 100
  expect_equal(im[3, 1], 100)

  # Check that attributes are preserved
  expect_equal(row_names(im), c(10L, 20L, 30L))
  expect_equal(col_names(im), c(100L, 200L, 300L, 400L))
})

test_that("getter and setter functions work", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m)

  # Set row names
  row_names(im) <- c(5L, 10L, 15L)
  expect_equal(row_names(im), c(5L, 10L, 15L))

  # Set column names
  col_names(im) <- c(50L, 100L, 150L, 200L)
  expect_equal(col_names(im), c(50L, 100L, 150L, 200L))

  # Validation in setters
  expect_error(
    row_names(im) <- c(1, 2, 3),
    "row_names must be integer or character vector"
  )

  expect_error(
    row_names(im) <- c(1L, 2L),
    "Length of row_names must match number of rows"
  )
})

test_that("character names work as dimnames", {
  m <- matrix(1:12, nrow = 3, ncol = 4)

  # Create with character names
  im <- matrix_int_names(m,
    row_names = c("row1", "row2", "row3"),
    col_names = c("col1", "col2", "col3", "col4")
  )

  expect_s3_class(im, "matrix_int_names")
  expect_equal(row_names(im), c("row1", "row2", "row3"))
  expect_equal(col_names(im), c("col1", "col2", "col3", "col4"))
  expect_true(is.character(row_names(im)))
  expect_true(is.character(col_names(im)))
})

test_that("subsetting by character names works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c("A", "B", "C"),
    col_names = c("W", "X", "Y", "Z")
  )

  # Subset by character row names
  result <- im[c("A", "C"), ]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(2, 4))
  expect_equal(row_names(result), c("A", "C"))
  expect_equal(as.vector(result[1, ]), c(1, 4, 7, 10))

  # Subset by character column names
  result <- im[, c("X", "Z")]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(3, 2))
  expect_equal(col_names(result), c("X", "Z"))

  # Subset by both
  result <- im[c("B", "C"), c("W", "Y")]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(2, 2))
  expect_equal(row_names(result), c("B", "C"))
  expect_equal(col_names(result), c("W", "Y"))
})

test_that("assignment with character names works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c("A", "B", "C"),
    col_names = c("W", "X", "Y", "Z")
  )

  # Assign using character names
  im["B", "Y"] <- 999
  expect_equal(im[2, 3], 999)

  # Check that attributes are preserved
  expect_equal(row_names(im), c("A", "B", "C"))
  expect_equal(col_names(im), c("W", "X", "Y", "Z"))
})

test_that("can mix integer and character names", {
  m <- matrix(1:12, nrow = 3, ncol = 4)

  # Integer row names, character column names
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c("W", "X", "Y", "Z")
  )

  expect_equal(row_names(im), c(10L, 20L, 30L))
  expect_equal(col_names(im), c("W", "X", "Y", "Z"))
  expect_true(is.integer(row_names(im)))
  expect_true(is.character(col_names(im)))

  # Can subset by both types
  result <- im[, , i_names = 10L, j_names = "X"]
  expect_equal(unname(result), 4)
})

test_that("row_names and col_names methods work with integer names", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c(100L, 200L, 300L, 400L)
  )

  # row_names() and col_names() should return integer names
  expect_equal(row_names(im), c(10L, 20L, 30L))
  expect_equal(col_names(im), c(100L, 200L, 300L, 400L))
  expect_true(is.integer(row_names(im)))
  expect_true(is.integer(col_names(im)))
})

test_that("row_names and col_names methods work with character names", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c("A", "B", "C"),
    col_names = c("W", "X", "Y", "Z")
  )

  # row_names() and col_names() should return character names
  expect_equal(row_names(im), c("A", "B", "C"))
  expect_equal(col_names(im), c("W", "X", "Y", "Z"))
  expect_true(is.character(row_names(im)))
  expect_true(is.character(col_names(im)))
})

test_that("row_names<- and col_names<- setters work", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m)

  # Set integer names using row_names<-
  row_names(im) <- c(5L, 10L, 15L)
  expect_equal(row_names(im), c(5L, 10L, 15L))
  expect_true(is.integer(row_names(im)))

  # Set character names using row_names<-
  row_names(im) <- c("A", "B", "C")
  expect_equal(row_names(im), c("A", "B", "C"))
  expect_true(is.character(row_names(im)))

  # Set integer names using col_names<-
  col_names(im) <- c(100L, 200L, 300L, 400L)
  expect_equal(col_names(im), c(100L, 200L, 300L, 400L))
  expect_true(is.integer(col_names(im)))

  # Set character names using col_names<-
  col_names(im) <- c("W", "X", "Y", "Z")
  expect_equal(col_names(im), c("W", "X", "Y", "Z"))
  expect_true(is.character(col_names(im)))
})

test_that("row_names and col_names setters clear names when set to NULL", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c(100L, 200L, 300L, 400L)
  )

  # Clear row names
  row_names(im) <- NULL
  expect_null(row_names(im))

  # Clear column names
  col_names(im) <- NULL
  expect_null(col_names(im))
})

test_that("row_names and col_names work with mixed name types", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c("W", "X", "Y", "Z")
  )

  # Integer row names, character column names
  expect_equal(row_names(im), c(10L, 20L, 30L))
  expect_equal(col_names(im), c("W", "X", "Y", "Z"))
  expect_true(is.integer(row_names(im)))
  expect_true(is.character(col_names(im)))
})

test_that("switching between integer and character names works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m, row_names = c(10L, 20L, 30L))

  # Start with integer names
  expect_equal(row_names(im), c(10L, 20L, 30L))
  expect_true(is.integer(row_names(im)))

  # Switch to character names
  row_names(im) <- c("A", "B", "C")
  expect_equal(row_names(im), c("A", "B", "C"))
  expect_true(is.character(row_names(im)))

  # Switch back to integer names
  row_names(im) <- c(100L, 200L, 300L)
  expect_equal(row_names(im), c(100L, 200L, 300L))
  expect_true(is.integer(row_names(im)))
})

test_that("logical subsetting works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c(100L, 200L, 300L, 400L)
  )

  # Logical row subsetting
  result <- im[c(TRUE, FALSE, TRUE), ]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(2, 4))
  expect_equal(row_names(result), c(10L, 30L))

  # Logical column subsetting
  result <- im[, c(TRUE, TRUE, FALSE, FALSE)]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(3, 2))
  expect_equal(col_names(result), c(100L, 200L))
})

test_that("negative indexing works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c(100L, 200L, 300L, 400L)
  )

  # Exclude rows
  result <- im[-2, ]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(2, 4))
  expect_equal(row_names(result), c(10L, 30L))

  # Exclude columns
  result <- im[, -c(1, 3)]
  expect_s3_class(result, "matrix_int_names")
  expect_equal(dim(result), c(3, 2))
  expect_equal(col_names(result), c(200L, 400L))
})

test_that("print method works", {
  m <- matrix(1:12, nrow = 3, ncol = 4)

  # With integer names
  im <- matrix_int_names(m,
    row_names = c(10L, 20L, 30L),
    col_names = c(100L, 200L, 300L, 400L)
  )
  expect_output(print(im), "Matrix with integer")
  expect_output(print(im), "3 x 4")
  expect_output(print(im), "10, 20, 30")
  expect_output(print(im), "100, 200, 300, 400")

  # With character names
  im2 <- matrix_int_names(m,
    row_names = c("A", "B", "C"),
    col_names = c("W", "X", "Y", "Z")
  )
  expect_output(print(im2), "Matrix with integer")
  expect_output(print(im2), "Character row names: A, B, C")
  expect_output(print(im2), "Character column names: W, X, Y, Z")
})

test_that("edge cases are handled", {
  # Empty matrix
  m <- matrix(numeric(0), nrow = 0, ncol = 0)
  im <- matrix_int_names(m, row_names = integer(0), col_names = integer(0))
  expect_s3_class(im, "matrix_int_names")
  expect_equal(dim(im), c(0, 0))

  # Single element matrix
  m <- matrix(42, nrow = 1, ncol = 1)
  im <- matrix_int_names(m, row_names = 5L, col_names = 10L)
  expect_s3_class(im, "matrix_int_names")
  expect_equal(im[, , i_names = 5L, j_names = 10L], 42)
})
