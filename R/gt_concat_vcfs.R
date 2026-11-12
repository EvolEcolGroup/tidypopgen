#' Concatenate VCF files of the same samples into a `gen_tibble`
#'
#' This function reads several *VCF* files with the same samples but different
#' loci (e.g. one file per chromosome, as distributed by many projects) into a
#' single [`gen_tibble`], as `bcftools concat` would, followed by
#' [gen_tibble()]. It needs no space on disk beyond the `gen_tibble`: each file
#' is read twice, once for its loci and samples, and once for its genotypes,
#' which are written straight into the backing file of the `gen_tibble`. To
#' combine files with different samples, use [gt_merge_vcfs()].
#'
#' As with [gen_tibble()], only biallelic SNPs are read from each file. The
#' loci of the `gen_tibble` are those of the first file, followed by those of
#' the second file, and so on (as in `bcftools concat`), so the files should
#' be given in the order of their chromosomes. The same chromosome can be
#' split across several files (e.g. in chunks), but duplicated loci (the same
#' chromosome and position) are an error, unless `allow_duplicates = TRUE`.
#'
#' All files must have the same samples, but they can be in a different order
#' in each file: the samples of the `gen_tibble` are in the order of the first
#' file, and the genotypes of the other files are reordered to match. Sample
#' names must be unique within each file.
#'
#' As in [gen_tibble()], the ploidy of each individual is estimated from the
#' first variants of each file; the ploidy of an individual in the
#' `gen_tibble` is its largest ploidy across files (e.g. with a file of the X
#' chromosome in which males are haploid, males are still diploid). Files in
#' which the genotypes of an individual are all missing in the first variants
#' are ignored for its ploidy (unless they all are).
#'
#' Unlike [gen_tibble()], by default (`names_as_int = TRUE`) the loci are
#' named with integers (their index in the `gen_tibble`), rather than with
#' their *VCF* IDs, to save memory with large datasets. With `names_as_int =
#' FALSE`, loci take their *VCF* IDs (a missing ID '.' gives a name made of
#' the chromosome and the position, as in [gen_tibble()]).
#'
#' @param files a character vector with the paths to two or more *VCF* files
#'   (.vcf or .vcf.gz, which are read faster if compressed with `bgzip`), in
#'   the order in which their loci should be in the `gen_tibble`.
#' @param backingfile the path, including the file name without extension, for
#'   the backing files of the `gen_tibble` (a .bk file for the genotypes, a
#'   .rds file for the FBM, and a .gt file for the `gen_tibble`). If NULL, the
#'   path of the first file is used, with "_concat" added to its name (and a
#'   version number if those files already exist, as in [gen_tibble()]).
#' @param valid_alleles a vector of valid allele values; it defaults to 'A','T',
#'   'C' and 'G'.
#' @param missing_alleles a vector of values in the *VCF* files that are
#'   considered as missing alleles ('.' is always missing).
#' @param allow_duplicates boolean: whether to allow duplicated loci (the same
#'   chromosome and position) or duplicated loci names in the `gen_tibble`
#'   (see [gen_tibble()]).
#' @param names_as_int boolean: whether the loci names are integers (the index
#'   of each locus in the `gen_tibble`) rather than the *VCF* IDs, which saves
#'   memory with many loci (see [gen_tibble()]). Note that it defaults to TRUE,
#'   unlike in [gen_tibble()] (where it defaults to FALSE): concatenating *VCF*
#'   files is mostly used for large datasets, in which the loci names take
#'   most of the memory of the `gen_tibble` and are rarely needed (loci are
#'   identified by their chromosome and position). Set it to FALSE to keep the
#'   *VCF* IDs. With integer names, the *VCF* IDs (with the chromosomes and
#'   positions) are written to a sidecar file next to the backing file, from
#'   which [gt_lookup_loci()] and [gt_restore_loci()] get them back.
#' @param chromosomes_as_int boolean: whether the chromosomes are integer
#'   codes rather than a factor, which saves memory with many chromosomes
#'   (e.g. RAD loci assembled de novo); the names are written to a sidecar
#'   file next to the backing file (see [gen_tibble()] and
#'   [gt_lookup_chromosomes()]). Default is FALSE.
#' @param n_cores the number of cores used to decompress bgzipped files and to
#'   parse the genotypes. Defaults to 1.
#' @param phased boolean: whether to store phased genotypes, giving a phased
#'   `gen_tibble` (see [gen_tibble_phased()]). All the files must then be
#'   phased: any unphased genotype (with '/') or individual that is not
#'   diploid is an error. Default is FALSE (the phase is ignored, as in
#'   [gen_tibble()]).
#' @param quiet boolean: whether to omit the report of the concatenation.
#' @returns a [`gen_tibble`] with the data of all files (saved, as with
#'   [gen_tibble()]).
#' @seealso [gt_merge_vcfs()] to merge *VCF* files with different samples.
#' @export
#' @examples
#' # two VCF files with the same individuals, one for each chromosome
#' example_gt <- load_example_gt("gen_tbl")
#' vcf_chr1 <- gt_as_vcf(select_loci(example_gt, 1:4),
#'   file = tempfile(fileext = ".vcf")
#' )
#' vcf_chr2 <- gt_as_vcf(select_loci(example_gt, 5:6),
#'   file = tempfile(fileext = ".vcf")
#' )
#' # concatenate them, keeping the VCF IDs as loci names (by default, unlike
#' # gen_tibble(), loci are named with integers to save memory)
#' concat_gt <- gt_concat_vcfs(c(vcf_chr1, vcf_chr2),
#'   names_as_int = FALSE,
#'   backingfile = tempfile("concat_")
#' )
#' concat_gt
#' show_loci(concat_gt)
gt_concat_vcfs <- function(
  files,
  backingfile = NULL,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  allow_duplicates = FALSE,
  names_as_int = TRUE,
  chromosomes_as_int = FALSE,
  n_cores = 1,
  phased = FALSE,
  quiet = FALSE
) {
  files <- check_vcf_files_args(files, valid_alleles, n_cores)
  n_cores <- as.integer(n_cores)
  backingfile <- default_vcfs_backingfile(backingfile, files, "_concat")
  n_files <- length(files)

  # first pass: the loci of each file (no genotypes), with its samples and
  # their ploidy
  # with integer names (or chromosomes), the original names are written to
  # sidecars next to the backing file (see R/loci_sidecar.R): the loci of
  # each file are added to the locus sidecar as the file is scanned (keyed
  # after those of the earlier files). The sidecars are removed if we stop
  sidecars <- vcfs_sidecar_paths(backingfile)
  sidecars_done <- FALSE
  on.exit(
    {
      unlink(sidecars$scan)
      if (!sidecars_done) unlink(c(sidecars$loci, sidecars$chromosomes))
    },
    add = TRUE
  )
  scans <- scan_vcf_files(files, names_as_int, n_cores, quiet,
    loci_file = if (names_as_int) sidecars$loci else ""
  )

  # the samples, in the order of the first file, and the row of each sample of
  # each file (an empty vector if the samples are in the same order)
  sample_id <- scans[[1]]$sample_names
  file_samples <- lapply(scans, function(x) x$sample_names)
  row_map <- concat_vcfs_row_maps(file_samples, files)
  # the rows (1-based) of the samples of each file
  file_rows <- lapply(row_map, function(x) {
    if (length(x) > 0) x + 1L else seq_along(sample_id)
  })
  # the ploidy of each individual: its largest ploidy across the files in
  # which it could be estimated. In a file where its genotypes were all
  # missing in the first variants, its ploidy was only guessed (from the
  # other individuals of that file), and that guess is only used if its
  # ploidy could not be estimated in any file (otherwise, e.g. a haploid
  # individual could get the diploid ploidy of the others)
  ploidy_known <- rep(NA_integer_, length(sample_id))
  ploidy_guess <- rep(NA_integer_, length(sample_id))
  for (f in seq_len(n_files)) {
    rows <- file_rows[[f]]
    guessed <- scans[[f]]$ploidy_guessed
    known <- ifelse(guessed, NA_integer_, scans[[f]]$ploidy)
    guess <- ifelse(guessed, scans[[f]]$ploidy, NA_integer_)
    ploidy_known[rows] <- pmax(ploidy_known[rows], known, na.rm = TRUE)
    ploidy_guess[rows] <- pmax(ploidy_guess[rows], guess, na.rm = TRUE)
  }
  ploidy <- ifelse(is.na(ploidy_known), ploidy_guess, ploidy_known)

  # the loci: those of each file, one file after the other. Each column is
  # taken from all the files, and then dropped from them, so that the loci
  # are not held twice in memory
  n_loci <- vapply(scans, function(x) length(x$loci_tbl$physical.pos), 1L)
  if (sum(n_loci) == 0) {
    stop("the VCF files have no biallelic SNPs")
  }
  # chromosomes get codes shared by all files, in order of first appearance
  chromosome_names <- vctrs::vec_unique(
    unlist(lapply(scans, function(x) x$chromosome_names))
  )
  # (the scans are held in an environment, so that their columns can be
  # dropped from them)
  scans_env <- new.env(parent = emptyenv())
  scans_env$scans <- scans
  rm(scans)
  take_loci_column <- function(column, fun) {
    x <- unlist(lapply(scans_env$scans, fun))
    for (f in seq_len(n_files)) {
      scans_env$scans[[f]]$loci_tbl[[column]] <- NULL
    }
    x
  }
  chromosome_index <- take_loci_column("chromosome", function(scan) {
    match(scan$chromosome_names, chromosome_names)[scan$loci_tbl$chromosome]
  })
  if (chromosomes_as_int) {
    # integer codes (see gen_tibble()), with their names in the sidecar
    chromosome_codes <- chromosome_int_codes_cpp(chromosome_names)
    loci_chromosome <- chromosome_codes[chromosome_index]
    write_chromosome_sidecar_cpp(
      sidecars$chromosomes, chromosome_codes,
      chromosome_names
    )
  } else {
    loci_chromosome <- chromosome_index_to_factor(
      chromosome_index,
      chromosome_names
    )
  }
  rm(chromosome_index)
  loci_position <- take_loci_column("physical.pos", function(scan) {
    scan$loci_tbl$physical.pos
  })
  if (names_as_int) {
    loci_name <- seq_along(loci_position)
  } else {
    # (missing IDs were filled with the chromosome and position when the
    # files were read)
    loci_name <- take_loci_column("marker.ID", function(scan) {
      scan$loci_tbl$marker.ID
    })
  }
  loci <- tibble::tibble(
    name = loci_name,
    chromosome = loci_chromosome,
    position = loci_position,
    genetic_dist = 0,
    allele_ref = take_loci_column("allele2", function(scan) {
      scan$loci_tbl$allele2
    }),
    allele_alt = take_loci_column("allele1", function(scan) {
      scan$loci_tbl$allele1
    })
  )
  rm(loci_name, loci_chromosome, loci_position)
  # a missing allele "." (e.g. the ALT of a monomorphic locus) is always
  # missing, as in gen_tibble()
  for (allele_col in c("allele_ref", "allele_alt")) {
    loci[[allele_col]][loci[[allele_col]] %in% "."] <- NA_character_
  }
  # validate the loci (before creating the FBM)
  loci <- validate_loci(loci,
    check_alphabet = TRUE,
    harmonise_loci = TRUE,
    check_duplicates = TRUE,
    allow_duplicates = allow_duplicates,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles,
    chromosomes_as_int = chromosomes_as_int
  )
  if (!quiet) {
    message(
      "concatenating ", n_files, " VCF files with ", length(sample_id),
      " samples: ", format(sum(n_loci), big.mark = ","),
      " loci (biallelic SNPs)"
    )
  }

  # second pass: the genotypes of each file, written straight into the
  # columns of its loci (all the genotypes of the gen_tibble are written)
  col_offset <- c(0L, cumsum(n_loci))
  x_gt <- write_vcfs_gen_tibble(
    files = files,
    loci = loci,
    sample_id = sample_id,
    ploidy = ploidy,
    row_offset = rep(0L, n_files),
    row_map = row_map,
    file_cols = function(f) {
      list(col = col_offset[f] + seq_len(n_loci[f]) - 1L, swap = logical(0))
    },
    file_ploidy = lapply(file_rows, function(x) ploidy[x]),
    fill_missing = FALSE,
    backingfile = backingfile,
    n_cores = n_cores,
    quiet = quiet,
    phased = phased,
    file_samples = file_samples
  )
  sidecars_done <- TRUE
  x_gt
}

# for each file, the 0-based row of each of its samples among the samples of
# the first file (or an empty vector if the samples are in the same order as
# in the first file); all files should have the same samples, with no
# duplicates
concat_vcfs_row_maps <- function(sample_names, files) {
  first_names <- sample_names[[1]]
  lapply(seq_along(sample_names), function(f) {
    these_names <- sample_names[[f]]
    if (anyDuplicated(these_names)) {
      stop("file ", f, " has duplicated sample names")
    }
    if (identical(these_names, first_names)) {
      return(integer(0))
    }
    absent <- setdiff(first_names, these_names)
    extra <- setdiff(these_names, first_names)
    if (length(absent) > 0 || length(extra) > 0) {
      problems <- c(
        if (length(absent) > 0) {
          paste0(
            "lacks samples of the first file (",
            paste(utils::head(absent, 5), collapse = ", "), ")"
          )
        },
        if (length(extra) > 0) {
          paste0(
            "has samples that are not in the first file (",
            paste(utils::head(extra, 5), collapse = ", "), ")"
          )
        }
      )
      stop(
        "all VCF files should have the same samples, but ", files[f], " ",
        paste(problems, collapse = " and ")
      )
    }
    match(these_names, first_names) - 1L
  })
}
