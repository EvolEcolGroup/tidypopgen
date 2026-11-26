test_that("TGENO files (made by EIGENSOFT) are read as their GENO files", {
  ext_dir <- system.file("extdata", package = "tidypopgen")
  tgeno_dir <- test_path("testdata", "eigensoft_tgeno")
  for (prefix in c("pop_a", "pop_b")) {
    geno_gt <- gen_tibble(file.path(ext_dir, paste0(prefix, ".geno")),
      backingfile = tempfile(), quiet = TRUE
    )
    # the TGENO file with the .snp and .ind files, with a .tgeno or a .geno
    # extension (the layout is read from the header), in lower or upper case
    tg_dir <- tempfile("tgeno_")
    dir.create(tg_dir)
    file.copy(file.path(ext_dir, paste0(prefix, c(".snp", ".ind"))), tg_dir)
    for (ext in c(".tgeno", ".geno", ".TGENO")) {
      tgeno_path <- file.path(tg_dir, paste0(prefix, ext))
      file.copy(file.path(tgeno_dir, paste0(prefix, ".tgeno")), tgeno_path)
      for (n_cores in c(1, 2)) {
        bk <- tempfile()
        tgeno_gt <- gen_tibble(tgeno_path,
          backingfile = bk, n_cores = n_cores, quiet = TRUE
        )
        expect_identical(show_genotypes(tgeno_gt), show_genotypes(geno_gt))
        expect_identical(show_loci(tgeno_gt), show_loci(geno_gt))
        expect_identical(tgeno_gt$id, geno_gt$id)
        expect_identical(tgeno_gt$population, geno_gt$population)
        expect_true(file.exists(paste0(bk, ".bk")))
      }
      unlink(tgeno_path)
    }
    # with no backingfile, the backing file is named after the .tgeno file
    file.copy(
      file.path(tgeno_dir, paste0(prefix, ".tgeno")),
      file.path(tg_dir, paste0(prefix, ".tgeno"))
    )
    tgeno_gt <- gen_tibble(file.path(tg_dir, paste0(prefix, ".tgeno")),
      quiet = TRUE
    )
    # (normalise both paths: on Windows the separators differ)
    expect_identical(
      normalizePath(gt_get_file_names(tgeno_gt)[2]),
      normalizePath(file.path(tg_dir, paste0(prefix, ".bk")))
    )
  }
})

test_that("GENO and TGENO files are read on several cores, in slabs", {
  set.seed(123)
  # n individuals (fewer than 192, so that GENO records are padded to 48
  # bytes, or more), and m loci, not multiples of 4 (the largest file is
  # decoded in several parallel tasks)
  for (dims in list(c(5, 1), c(1, 9), c(37, 1037), c(203, 13001))) {
    n <- dims[1]
    m <- dims[2]
    # packed values: the number of copies of the first allele, or 3 (missing)
    values <- matrix(
      sample(0:3, n * m, replace = TRUE, prob = c(0.3, 0.3, 0.3, 0.1)),
      nrow = n, ncol = m
    )
    # the dosages of the second allele (as doubles, as returned by the FBM)
    dosages <- 2 - values
    dosages[values == 3] <- NA
    prefix <- tempfile("packed_")
    write_packed_snp_ind(prefix, n = n, m = m)
    geno_file <- write_packed_geno(values, paste0(prefix, ".geno"))
    tgeno_file <- write_packed_geno(values, paste0(prefix, ".tgeno"),
      transposed = TRUE
    )
    for (n_cores in c(1, 2)) {
      for (file in c(geno_file, tgeno_file)) {
        gt <- gen_tibble(file,
          backingfile = tempfile(), n_cores = n_cores, quiet = TRUE
        )
        expect_identical(show_genotypes(gt), dosages)
      }
      # a TGENO file read in slabs of loci of 64 or 128 bytes (the buffer is
      # much smaller than the file), or of nearly half the file (the buffer
      # is just smaller than the file), or in one go (it just fits)
      record_bytes <- ceiling(m / 4)
      for (max_buffer_bytes in c(1, 2 * 128 * n, n * record_bytes - 1,
                                 n * record_bytes)) {
        fbm <- fbm_read_packedancestry(tgeno_file,
          n_indiv = n, n_snp = m, backingfile = tempfile(),
          n_cores = n_cores, max_buffer_bytes = max_buffer_bytes
        )
        expect_identical(fbm[, , drop = FALSE], dosages)
        # (and the FBM saved in its .rds file)
        fbm_rds <- readRDS(bigstatsr::sub_bk(fbm$backingfile, ".rds"))
        expect_identical(fbm_rds[, , drop = FALSE], dosages)
      }
    }
  }
})

test_that("bad packedancestry files are caught", {
  ext_dir <- system.file("extdata", package = "tidypopgen")
  tgeno_dir <- test_path("testdata", "eigensoft_tgeno")
  bad_dir <- tempfile("bad_geno_")
  dir.create(bad_dir)
  file.copy(file.path(ext_dir, c("pop_a.snp", "pop_a.ind")), bad_dir)
  bad_path <- file.path(bad_dir, "pop_a.geno")
  # each bad file gives an error, and leaves no backing file behind
  expect_bad_geno <- function(bytes, message) {
    writeBin(bytes, bad_path)
    bk <- tempfile()
    expect_error(
      gen_tibble(bad_path, backingfile = bk, quiet = TRUE),
      message
    )
    expect_false(file.exists(paste0(bk, ".bk")))
  }
  geno_bytes <- readBin(file.path(ext_dir, "pop_a.geno"), "raw", n = 1e5)
  tgeno_bytes <- readBin(file.path(tgeno_dir, "pop_a.tgeno"), "raw", n = 1e5)
  for (bytes in list(geno_bytes, tgeno_bytes)) {
    # truncated (a byte missing, or only the header)
    expect_bad_geno(bytes[-length(bytes)], "truncated or corrupted")
    expect_bad_geno(bytes[1:48], "truncated or corrupted")
  }
  # numbers of individuals and of loci that do not match the .ind and .snp
  header_with <- function(header) {
    header <- charToRaw(header)
    c(header, raw(48 - length(header)), tgeno_bytes[-(1:48)])
  }
  expect_bad_geno(
    header_with("TGENO       6      16 0 0"),
    "Number of individuals in geno file does not match"
  )
  expect_bad_geno(
    header_with("TGENO       5      17 0 0"),
    "Number of variants in geno file does not match"
  )
  expect_bad_geno(header_with("TGENO 5"), "malformed")
  # not a packedancestrymap file
  expect_bad_geno(c(as.raw(0), tgeno_bytes[-1]), "not a packedancestrymap")
  expect_bad_geno(raw(0), "not a packedancestrymap")
  # an EIGENSTRAT (text) .geno file
  expect_bad_geno(
    charToRaw(strrep("00102\n", 16)),
    "looks like an EIGENSTRAT \\(text\\) .geno file"
  )
  # extra bytes at the end are ignored, with a warning
  for (bytes in list(geno_bytes, tgeno_bytes)) {
    writeBin(c(bytes, as.raw(c(0, 255))), bad_path)
    expect_warning(
      extra_gt <- gen_tibble(bad_path, backingfile = tempfile(), quiet = TRUE),
      "the extra bytes at its end were ignored"
    )
    pop_a_gt <- gen_tibble(file.path(ext_dir, "pop_a.geno"),
      backingfile = tempfile(), quiet = TRUE
    )
    expect_identical(show_genotypes(extra_gt), show_genotypes(pop_a_gt))
  }
})

test_that("packedancestry files are decoded with any code table", {
  # a random code table, which is not made of codes of 2-bit genotypes, so
  # that TGENO files are decoded one byte at a time (as on big-endian
  # platforms), checked against the bytes the table gives
  set.seed(42)
  n <- 21
  m <- 203
  values <- matrix(sample(0:3, n * m, replace = TRUE), nrow = n, ncol = m)
  prefix <- tempfile("packed_")
  geno_file <- write_packed_geno(values, paste0(prefix, ".geno"))
  tgeno_file <- write_packed_geno(values, paste0(prefix, ".tgeno"),
    transposed = TRUE
  )
  # TGENO: byte j %/% 4 of the record of individual i holds the genotype of
  # locus j (counting from 0), in position j %% 4 (from the highest bits);
  # GENO: byte i %/% 4 of the record of locus j holds that of individual i,
  # in position i %% 4
  tgeno_records <- matrix(
    as.integer(readBin(tgeno_file, "raw", n = 1e6)[-(1:48)]),
    nrow = ceiling(m / 4)
  )
  geno_record_bytes <- max(48, ceiling(n / 4))
  geno_records <- matrix(
    as.integer(readBin(geno_file, "raw", n = 1e6)[-seq_len(geno_record_bytes)]),
    nrow = geno_record_bytes
  )
  for (tab in list(
    matrix(as.raw(sample(0:255, 4 * 256, replace = TRUE)), nrow = 4),
    get_packedancestry_code()
  )) {
    expected_tgeno <- matrix(raw(n * m), nrow = n, ncol = m)
    for (j in seq_len(m)) {
      expected_tgeno[, j] <- tab[
        (j - 1) %% 4 + 1,
        tgeno_records[(j - 1) %/% 4 + 1, ] + 1
      ]
    }
    expected_geno <- matrix(raw(n * m), nrow = n, ncol = m)
    for (i in seq_len(n)) {
      expected_geno[i, ] <- tab[
        (i - 1) %% 4 + 1,
        geno_records[(i - 1) %/% 4 + 1, ] + 1
      ]
    }
    for (n_threads in c(1, 2)) {
      for (transposed in c(FALSE, TRUE)) {
        for (max_buffer_bytes in c(1, 2^31)) {
          bk <- tempfile()
          fbm <- bigstatsr::FBM.code256(n, m,
            code = bigsnpr::CODE_012, backingfile = bk
          )
          read_packedancestry_to_fbm(
            if (transposed) tgeno_file else geno_file, fbm, tab,
            transposed = transposed, n_threads = n_threads,
            max_buffer_bytes = max_buffer_bytes
          )
          expect_identical(
            readBin(paste0(bk, ".bk"), "raw", n = n * m),
            as.vector(if (transposed) expected_tgeno else expected_geno)
          )
        }
      }
    }
  }
})
