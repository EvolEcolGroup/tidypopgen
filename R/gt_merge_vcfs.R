#' Merge several VCF files into a `gen_tibble`
#'
#' This function reads several *VCF* files with different samples (e.g. one
#' per population, or per sequencing batch) into a single [`gen_tibble`], as
#' `bcftools merge` would, followed by [gen_tibble()]. It is much faster than
#' reading each file with [gen_tibble()] and merging them with
#' [rbind()][rbind.gen_tbl] (one pair at a time), and needs no space on disk
#' beyond the merged `gen_tibble`: each file is read twice, once for its loci
#' (to match them across files) and once for its genotypes, which are written
#' straight into the merged backing file.
#'
#' As with [gen_tibble()], only biallelic SNPs are read from each file. Loci
#' are matched across files by chromosome and position (as in `bcftools
#' merge`), so chromosomes need to be named in the same way in all files (a mix
#' of e.g. 'chr1' and '1' would match no loci). A locus is merged with a site
#' at the same chromosome and position if their alleles agree: the same REF
#' and ALT alleles, or a missing allele (e.g. the ALT '.' of a monomorphic
#' locus), which takes the allele of the other files. Loci whose alleles do not
#' agree are handled as follows:
#' - loci with the same REF but a different ALT allele (e.g. A/C and A/G):
#'   `bcftools merge` would make a multiallelic site, which a `gen_tibble` can
#'   not store, so by default the site is dropped (`alt_mismatch = "drop"`);
#'   with `alt_mismatch = "split"` (as `bcftools merge -m none`), the loci are
#'   kept as separate loci at the same position (which needs `allow_duplicates
#'   = TRUE`), each with missing genotypes for the samples of the other files.
#' - loci with a different REF allele (e.g. A/C and C/A): the REF allele should
#'   be the base of the reference genome at that position, so `bcftools merge`
#'   stops with an error, as does this function by default (`ref_mismatch =
#'   "error"`, listing the first few mismatches). With `ref_mismatch = "drop"`
#'   such sites are dropped; with `ref_mismatch = "swap"`, loci whose alleles
#'   are swapped (A/C and C/A) are merged, swapping the genotypes of the later
#'   file (as [rbind()][rbind.gen_tbl] does), and other mismatches are handled
#'   as with `alt_mismatch`. Swapping is needed e.g. for *VCF* files exported
#'   by PLINK 1.9, which by default set the REF allele to the major allele.
#'
#' Strands are not flipped (all *VCF* files are on the forward strand). Loci
#' at the same position within a file are never merged with each other (as in
#' `bcftools merge`).
#'
#' The merged `gen_tibble` has the samples of the first file, followed by those
#' of the second file, and so on. By default (`variants = "intersection"`),
#' it keeps the sites found in all files, as [rbind()][rbind.gen_tbl] does;
#' with `variants = "union"`, it keeps all sites (as `bcftools merge`), and the
#' samples of the
#' files without a site get missing genotypes there, or homozygous reference
#' genotypes with `missing_to_ref = TRUE` (as `bcftools merge -0`). Sites are
#' ordered by chromosome (in order of first appearance across the files), then
#' by position. Each locus takes its name (*VCF* ID) from the first file where
#' it is found (a missing ID '.' gives a name made of the chromosome and the
#' position, as in [gen_tibble()]).
#'
#' Sample names must be unique across files: if they are not, the function
#' stops, unless `force_samples = TRUE`, in which case the duplicated names are
#' prefixed by the index of their file, as `bcftools merge --force-samples`
#' does (e.g. sample 'S3' of the second file becomes '2:S3').
#'
#' As in [gen_tibble()], the ploidy of each individual is estimated from the
#' first variants of its file. The individuals of different files can have
#' different ploidies.
#'
#' @param files a character vector with the paths to two or more *VCF* files
#'   (.vcf or .vcf.gz, which are read faster if compressed with `bgzip`).
#' @param variants which sites to keep: "intersection" (the default) keeps the
#'   sites found in all files; "union" keeps all sites (see details).
#' @param missing_to_ref boolean: with `variants = "union"`, whether the
#'   genotypes of the samples of a file at sites that are not in that file are
#'   homozygous for the reference allele (as with `bcftools merge -0`), rather
#'   than missing (the default).
#' @param alt_mismatch what to do with loci at the same position that have the
#'   same REF but a different ALT allele: "drop" (the default) drops the site,
#'   "split" keeps them as separate loci (see details).
#' @param ref_mismatch what to do with loci at the same position that have a
#'   different REF allele: "error" (the default) stops with an error, "drop"
#'   drops the site, "swap" merges loci whose alleles are swapped, swapping
#'   their genotypes (see details).
#' @param force_samples boolean: whether to allow the same sample name in
#'   different files, prefixing the duplicated names by the index of their file
#'   (see details). Defaults to FALSE.
#' @param backingfile the path, including the file name without extension, for
#'   the backing files of the merged `gen_tibble` (a .bk file for the
#'   genotypes, a .rds file for the FBM, and a .gt file for the `gen_tibble`).
#'   If NULL, the path of the first file is used, with "_merged" added to its
#'   name (and a version number if those files already exist, as in
#'   [gen_tibble()]).
#' @param valid_alleles a vector of valid allele values; it defaults to 'A','T',
#'   'C' and 'G'.
#' @param missing_alleles a vector of values in the *VCF* files that are
#'   considered as missing alleles ('.' is always missing).
#' @param allow_duplicates boolean: whether to allow duplicated loci (the same
#'   chromosome and position) or duplicated loci names in the merged
#'   `gen_tibble` (see [gen_tibble()]).
#' @param names_as_int boolean: whether the loci names are integers (the index
#'   of each locus in the merged `gen_tibble`) rather than the *VCF* IDs, which
#'   saves memory with many loci (see [gen_tibble()]).
#' @param n_cores the number of cores used to decompress bgzipped files and to
#'   parse the genotypes. Defaults to 1.
#' @param quiet boolean: whether to omit the report of the merge.
#' @returns a [`gen_tibble`] with the merged data (saved, as with
#'   [gen_tibble()]).
#' @seealso [rbind()][rbind.gen_tbl] to merge two `gen_tibble` objects (e.g.
#'   read from files of other formats).
#' @export
#' @examples
#' vcf_a <- system.file("extdata", "pop_a.vcf", package = "tidypopgen")
#' vcf_b <- system.file("extdata", "pop_b.vcf", package = "tidypopgen")
#' # the REF and ALT alleles of some loci are swapped between these two files
#' # (which a VCF from PLINK can have), so we swap them back to merge them
#' merged_gt <- gt_merge_vcfs(c(vcf_a, vcf_b),
#'   ref_mismatch = "swap",
#'   backingfile = tempfile("merged_")
#' )
#' merged_gt
#' show_loci(merged_gt)
gt_merge_vcfs <- function(
  files,
  variants = c("intersection", "union"),
  missing_to_ref = FALSE,
  alt_mismatch = c("drop", "split"),
  ref_mismatch = c("error", "drop", "swap"),
  force_samples = FALSE,
  backingfile = NULL,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  allow_duplicates = FALSE,
  names_as_int = FALSE,
  n_cores = 1,
  quiet = FALSE
) {
  variants <- match.arg(variants)
  alt_mismatch <- match.arg(alt_mismatch)
  ref_mismatch <- match.arg(ref_mismatch)
  if (!is.character(files) || length(files) < 2) {
    stop(
      "'files' should be the paths to two or more VCF files ",
      "(use gen_tibble() to read a single file)"
    )
  }
  if (anyNA(files) || any(!file.exists(files))) {
    stop(
      "these VCF files do not exist: ",
      paste(files[is.na(files) | !file.exists(files)], collapse = ", ")
    )
  }
  if (!all(grepl("\\.(vcf|vcf\\.gz)$", files, ignore.case = TRUE))) {
    stop("'files' should be VCF files (with extension .vcf or .vcf.gz)")
  }
  files <- normalizePath(files)
  if (anyDuplicated(files)) {
    stop("'files' contains the same file more than once")
  }
  # check that valid alleles does not contain zero
  if ("0" %in% valid_alleles) {
    stop(paste(
      "'0' cannot be a valid allele",
      "(it is the default missing allele value!)"
    ))
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || is.na(n_cores) ||
    n_cores < 1) {
    stop("'n_cores' should be a positive integer")
  }
  n_cores <- as.integer(n_cores)
  if (missing_to_ref && variants == "intersection" && !quiet) {
    message(
      "'missing_to_ref' is only used with variants = 'union' (with ",
      "'intersection', all files have all sites)"
    )
  }
  # sort out the backing file (by default, next to the first file)
  if (is.null(backingfile)) {
    backingfile <- paste0(
      sub("\\.(vcf\\.gz|vcf)$", "", files[1], ignore.case = TRUE),
      "_merged"
    )
  }
  backingfile <- change_duplicated_file_name(backingfile)
  n_files <- length(files)

  # first pass: the loci of each file (no genotypes), with its samples and
  # their ploidy (locus names are not read if we use integers)
  scans <- vector("list", n_files)
  for (f in seq_len(n_files)) {
    if (!quiet) {
      message("reading the loci of ", files[f])
    }
    scans[[f]] <- vcf_to_bk(files[f],
      bk_file = "", read_names = !names_as_int,
      n_threads = n_cores, read_genotypes = FALSE
    )
    if (length(scans[[f]]$sample_names) == 0) {
      stop("the VCF ", files[f], " has no samples (or no variants)")
    }
  }

  # the samples (in the order of the files)
  sample_id <- merge_vcfs_sample_names(
    lapply(scans, function(x) x$sample_names),
    force_samples = force_samples
  )
  n_indiv <- vapply(scans, function(x) length(x$sample_names), integer(1))
  ploidy <- unlist(lapply(scans, function(x) x$ploidy))
  if (anyNA(ploidy) || any(ploidy <= 0)) {
    stop("the ploidy of each individual should be a positive integer")
  }
  max_ploidy <- max(ploidy)

  # chromosomes: each file numbers its own (unique) chromosome names, so we
  # give them codes shared by all files, in order of first appearance across
  # the files
  chromosome_names <- vctrs::vec_unique(
    unlist(lapply(scans, function(x) x$chromosome_names))
  )
  chromosome <- lapply(scans, function(x) {
    match(x$chromosome_names, chromosome_names)[x$loci_tbl$chromosome]
  })
  # alleles, with all missing alleles as NA (a missing allele "." is always
  # missing, as in gen_tibble())
  all_missing <- union(c(missing_alleles, "0"), ".")
  allele_as_na <- function(x) {
    x[x %in% all_missing] <- NA_character_
    x
  }

  # match the loci of the files into sites
  sites <- vcf_merge_sites_cpp(
    chromosome = chromosome,
    position = lapply(scans, function(x) x$loci_tbl$physical.pos),
    allele_ref = lapply(scans, function(x) allele_as_na(x$loci_tbl$allele2)),
    allele_alt = lapply(scans, function(x) allele_as_na(x$loci_tbl$allele1)),
    split_alt = alt_mismatch == "split",
    ref_rule = match(ref_mismatch, c("error", "drop", "swap")) - 1L
  )
  if (sites$ref_mismatch$n > 0) {
    stop_ref_mismatch(sites, scans, chromosome_names, files)
  }

  # the sites we keep, ordered by chromosome (in order of first appearance)
  # and position (sites at the same position stay in the order in which they
  # were found)
  keep <- sites$status == 0L
  if (variants == "intersection") {
    keep <- keep & sites$n_files == n_files
  }
  keep <- which(keep)
  keep <- keep[order(
    sites$chromosome[keep], sites$position[keep],
    method = "radix"
  )]
  if (length(keep) == 0) {
    stop(
      "there are no loci in common between the VCF files ",
      "(check that chromosomes are named in the same way in all files)"
    )
  }
  if (!quiet) {
    report_vcf_merge(sites, keep, variants, n_files, sum(n_indiv))
  }

  # the loci table of the merged gen_tibble; its chromosomes are only those
  # with merged loci (as in gen_tibble(), whose chromosomes are those with
  # biallelic SNPs), renumbered in the same order
  chromosome_used <- rep(FALSE, length(chromosome_names))
  chromosome_used[sites$chromosome[keep]] <- TRUE
  chromosome_new_index <- cumsum(chromosome_used)
  loci_chromosome <- chromosome_index_to_factor(
    chromosome_new_index[sites$chromosome[keep]],
    chromosome_names[chromosome_used]
  )
  loci_position <- sites$position[keep]
  if (names_as_int) {
    loci_name <- seq_along(keep)
  } else {
    # each locus takes the name it has in the first file where it is found
    first_file <- sites$first_file[keep]
    first_locus <- sites$first_locus[keep]
    loci_name <- character(length(keep))
    for (f in unique(first_file)) {
      in_f <- which(first_file == f)
      loci_name[in_f] <- scans[[f]]$loci_tbl$marker.ID[first_locus[in_f]]
    }
    loci_name <- fill_missing_loci_names(
      loci_name, loci_chromosome,
      loci_position
    )
  }
  loci <- tibble::tibble(
    name = loci_name,
    chromosome = loci_chromosome,
    position = loci_position,
    genetic_dist = 0,
    allele_ref = sites$allele_ref[keep],
    allele_alt = sites$allele_alt[keep]
  )
  # validate the loci (before creating the FBM)
  loci <- validate_loci(loci,
    check_alphabet = TRUE,
    harmonise_loci = TRUE,
    check_duplicates = TRUE,
    allow_duplicates = allow_duplicates,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles
  )
  indiv_meta <- validate_indiv_meta(data.frame(id = sample_id))

  # the column of the merged FBM (0-based) of each site (-1 if not kept)
  site_col <- rep(-1L, length(sites$status))
  site_col[keep] <- seq_along(keep) - 1L

  # the merged FBM: missing genotypes are its largest ploidy + 1. With
  # variants = "union", the samples of a file have missing genotypes at the
  # sites that are not in that file (or 0, with missing_to_ref = TRUE, which
  # is how an FBM starts), so we start with all genotypes missing; with
  # "intersection", all genotypes are written in the second pass
  missing_value <- max_ploidy + 1
  code256 <- rep(NA_real_, 256)
  code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)
  fill_missing <- variants == "union" && !missing_to_ref
  file_backed_matrix <- bigstatsr::FBM.code256(
    nrow = sum(n_indiv),
    ncol = length(keep),
    code = code256,
    backingfile = backingfile,
    init = if (fill_missing) missing_value else NULL
  )
  bk_file <- file_backed_matrix$backingfile
  # remove the backing file if we stop before creating the gen_tibble
  fbm_done <- FALSE
  on.exit(
    if (!fbm_done) {
      rm(file_backed_matrix)
      remove_fbm_files(c(bk_file, bigstatsr::sub_bk(bk_file, ".rds")))
    },
    add = TRUE
  )

  # second pass: the genotypes of each file, written straight into the
  # columns of its loci, in the rows of its samples
  row_offset <- c(0L, cumsum(n_indiv))
  for (f in seq_len(n_files)) {
    if (!quiet) {
      message("reading the genotypes of ", files[f])
    }
    vcf_merge_genotypes(files[f],
      BM = file_backed_matrix,
      row_offset = row_offset[f],
      col = site_col[sites$site[[f]]],
      swap = sites$swap[[f]],
      ploidy = scans[[f]]$ploidy,
      missing_value = missing_value,
      n_threads = n_cores
    )
  }
  # save the FBM (before its md5 sum is taken by new_vctrs_bigsnp())
  fbm_path <- bigstatsr::sub_bk(bk_file, ".rds")
  saveRDS(file_backed_matrix, fbm_path)

  indiv_meta$genotypes <- new_vctrs_bigsnp(
    fbm_obj = file_backed_matrix,
    fbm_file = fbm_path,
    loci = loci,
    indiv_id = indiv_meta$id,
    ploidy = max_ploidy,
    fbm_ploidy = ploidy
  )
  x_gt <- tibble::new_tibble(indiv_meta, class = "gen_tbl")
  fbm_done <- TRUE
  rm(file_backed_matrix)

  # check alleles (as gen_tibble() does)
  if (has_genotypes_without_alleles(x_gt)) {
    # (dropping the gen_tibble first, so that its backing file can be
    # unmapped, see remove_fbm_files())
    gt_files <- gt_get_file_names(x_gt)
    rm(x_gt)
    remove_fbm_files(gt_files)
    stop_genotypes_without_alleles()
  }
  gt_save_light(x_gt, quiet = quiet) # nolint
  x_gt
}

# the sample names of the merged VCFs, from a list of the sample names of each
# file: duplicated names are an error, unless force_samples is TRUE, in which
# case names already used by an earlier file are prefixed by the index of
# their file (as bcftools merge --force-samples does, e.g. "2:S3")
merge_vcfs_sample_names <- function(sample_names, force_samples = FALSE) {
  # names duplicated within a file are always an error
  for (f in seq_along(sample_names)) {
    if (anyDuplicated(sample_names[[f]])) {
      stop("file ", f, " has duplicated sample names")
    }
  }
  all_names <- unlist(sample_names)
  if (!anyDuplicated(all_names)) {
    return(all_names)
  }
  if (!force_samples) {
    stop(
      "the same sample names are found in more than one file: ",
      paste(utils::head(unique(all_names[duplicated(all_names)]), 5),
        collapse = ", "
      ),
      ". Use force_samples = TRUE to rename them."
    )
  }
  seen <- character(0)
  for (f in seq_along(sample_names)) {
    dup <- sample_names[[f]] %in% seen
    sample_names[[f]][dup] <- paste0(f, ":", sample_names[[f]][dup])
    seen <- c(seen, sample_names[[f]])
  }
  all_names <- unlist(sample_names)
  if (anyDuplicated(all_names)) {
    stop(
      "renaming the duplicated sample names created more duplicated names: ",
      paste(utils::head(unique(all_names[duplicated(all_names)]), 5),
        collapse = ", "
      )
    )
  }
  all_names
}

# stop because the REF alleles of loci at the same position differ between
# files, listing the first few
stop_ref_mismatch <- function(sites, scans, chromosome_names, files) {
  mm <- sites$ref_mismatch
  describe <- function(f, k) {
    paste0(
      scans[[f]]$loci_tbl$allele2[k], "/",
      scans[[f]]$loci_tbl$allele1[k], " in ", basename(files[f])
    )
  }
  examples <- vapply(seq_along(mm$site), function(i) {
    s <- mm$site[i]
    paste0(
      chromosome_names[sites$chromosome[s]], ":", sites$position[s], " (",
      describe(sites$first_file[s], sites$first_locus[s]), ", ",
      describe(mm$file[i], mm$locus[i]), ")"
    )
  }, character(1))
  stop(
    "the REF alleles of ", format(mm$n, big.mark = ","), " loci differ ",
    "between files (REF/ALT shown), e.g. ",
    paste(examples, collapse = "; "),
    ". Use ref_mismatch = 'swap' to merge loci with swapped alleles, or ",
    "ref_mismatch = 'drop' to drop these loci.",
    call. = FALSE
  )
}

# report the numbers of sites of a merge
report_vcf_merge <- function(sites, keep, variants, n_files, n_indiv) {
  n_sites <- length(sites$status)
  fmt <- function(x) format(x, big.mark = ",")
  message(
    "merging ", n_files, " VCF files with ", fmt(n_indiv), " samples:\n",
    "  ", fmt(n_sites), " sites (biallelic SNPs), of which ",
    fmt(sum(sites$n_files == n_files)), " are in all files\n",
    "  dropped: ", fmt(sum(sites$status == 1L)),
    " with alleles that would make a multiallelic site, ",
    fmt(sum(sites$status == 2L)), " with different REF alleles\n",
    "  loci with swapped alleles: ",
    fmt(sum(vapply(sites$swap, sum, numeric(1)))), "\n",
    "  kept: ", fmt(length(keep)), " loci (", variants, ")"
  )
}
