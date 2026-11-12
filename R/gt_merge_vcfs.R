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
#' locus), which takes the allele of the other files (a site with the same
#' alleles is preferred to one with a missing allele). Loci whose alleles do
#' not agree are handled as follows:
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
#' `bcftools merge`); if a file has several loci at the same position (e.g. a
#' multiallelic site split into biallelic records, as by `bcftools norm -m-`),
#' each is merged with the site of the other files that has its alleles, if
#' any, before the others are handled as above.
#'
#' The merged `gen_tibble` has the samples of the first file, followed by those
#' of the second file, and so on. By default (`variants = "intersection"`),
#' it keeps the sites found in all files, as [rbind()][rbind.gen_tbl] does;
#' with `variants = "union"`, it keeps all sites (as `bcftools merge`), and the
#' samples of the
#' files without a site get missing genotypes there, or homozygous reference
#' genotypes with `missing_to_ref = TRUE` (as `bcftools merge -0`). Sites are
#' ordered by chromosome (in order of first appearance across the files), then
#' by position.
#'
#' Unlike [gen_tibble()], by default (`names_as_int = TRUE`) the loci are
#' named with integers (their index in the merged `gen_tibble`), rather than
#' with their *VCF* IDs, to save memory with large datasets. With
#' `names_as_int = FALSE`, each locus takes its name (*VCF* ID) from the first
#' file where it is found (a missing ID '.' gives a name made of the
#' chromosome and the position, as in [gen_tibble()]).
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
#'   saves memory with many loci (see [gen_tibble()]). Note that it defaults
#'   to TRUE, unlike in [gen_tibble()] (where it defaults to FALSE): merging
#'   *VCF* files is mostly used for large datasets, in which the loci names
#'   take most of the memory of the `gen_tibble` and are rarely needed (loci
#'   are identified by their chromosome and position). Set it to FALSE to keep
#'   the *VCF* IDs. With integer names, the *VCF* IDs (with the chromosomes
#'   and positions) are written to a sidecar file next to the backing file,
#'   from which [gt_lookup_loci()] and [gt_restore_loci()] get them back.
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
#'   [gen_tibble()]). The alleles of phased genotypes are never swapped (as in
#'   [rbind()][rbind.gen_tbl]), so `ref_mismatch = "swap"` is not available.
#' @param quiet boolean: whether to omit the report of the merge.
#' @returns a [`gen_tibble`] with the merged data (saved, as with
#'   [gen_tibble()]).
#' @seealso [gt_concat_vcfs()] to combine *VCF* files with the same samples
#'   and different loci (e.g. one file per chromosome);
#'   [rbind()][rbind.gen_tbl] to merge two `gen_tibble` objects (e.g. read
#'   from files of other formats).
#' @export
#' @examples
#' vcf_a <- system.file("extdata", "pop_a.vcf", package = "tidypopgen")
#' vcf_b <- system.file("extdata", "pop_b.vcf", package = "tidypopgen")
#' # the REF and ALT alleles of some loci are swapped between these two files
#' # (which a VCF from PLINK can have), so we swap them back to merge them;
#' # we also keep the VCF IDs as loci names (by default, unlike gen_tibble(),
#' # loci are named with integers to save memory)
#' merged_gt <- gt_merge_vcfs(c(vcf_a, vcf_b),
#'   ref_mismatch = "swap",
#'   names_as_int = FALSE,
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
  names_as_int = TRUE,
  chromosomes_as_int = FALSE,
  n_cores = 1,
  phased = FALSE,
  quiet = FALSE
) {
  variants <- match.arg(variants)
  alt_mismatch <- match.arg(alt_mismatch)
  ref_mismatch <- match.arg(ref_mismatch)
  files <- check_vcf_files_args(files, valid_alleles, n_cores)
  if (isTRUE(phased) && ref_mismatch == "swap") {
    stop(
      "the alleles of phased genotypes are never swapped (as in rbind()): ",
      "use ref_mismatch = 'error' or 'drop'"
    )
  }
  n_cores <- as.integer(n_cores)
  if (missing_to_ref && variants == "intersection" && !quiet) {
    message(
      "'missing_to_ref' is only used with variants = 'union' (with ",
      "'intersection', all files have all sites)"
    )
  }
  backingfile <- default_vcfs_backingfile(backingfile, files, "_merged")
  n_files <- length(files)

  # with integer names (or chromosomes), the original names are written to
  # sidecars next to the backing file (see R/loci_sidecar.R). Each scan
  # writes the loci of its file to a temporary sidecar, from which the loci
  # that create new sites are appended to a sidecar of the sites (keyed by
  # site), so that only the loci of one file and the sites are on disk at a
  # time; the sites that are kept are gathered at the end, in the order of
  # the merged loci. The sidecars are removed if we stop
  sidecars <- vcfs_sidecar_paths(backingfile)
  sidecars_done <- FALSE
  on.exit(
    {
      unlink(c(sidecars$scan, sidecars$sites))
      if (!sidecars_done) unlink(c(sidecars$loci, sidecars$chromosomes))
    },
    add = TRUE
  )
  n_sites <- 0

  # first pass: the loci of each file (no genotypes), with its samples and
  # their ploidy, matched into sites one file at a time (so that only the
  # loci of one file are in memory at a time). Missing alleles (a missing
  # allele "." is always missing, as in gen_tibble()) can be any allele
  merged_sites <- merge_sites_new(
    split_alt = alt_mismatch == "split",
    ref_rule = match(ref_mismatch, c("error", "drop", "swap")) - 1L,
    missing_alleles = union(c(missing_alleles, "0"), ".")
  )
  sample_names <- vector("list", n_files)
  file_ploidy <- vector("list", n_files)
  # the names of the sites (the names of the loci that created them)
  site_names <- vector("list", n_files)
  # chromosomes get codes shared by all files, in order of first appearance
  chromosome_names <- character(0)
  for (f in seq_len(n_files)) {
    scan <- scan_vcf_file(files[f], names_as_int, n_cores, quiet,
      loci_file = if (names_as_int) sidecars$scan else ""
    )
    sample_names[[f]] <- scan$sample_names
    file_ploidy[[f]] <- scan$ploidy
    # check the samples as soon as possible
    merge_vcfs_sample_names(sample_names[seq_len(f)], force_samples)
    chromosome_names <- vctrs::vec_unique(
      c(chromosome_names, scan$chromosome_names)
    )
    new_loci <- merge_sites_add_file(merged_sites,
      chromosome = match(scan$chromosome_names, chromosome_names)[
        scan$loci_tbl$chromosome
      ],
      position = scan$loci_tbl$physical.pos,
      allele_ref = scan$loci_tbl$allele2,
      allele_alt = scan$loci_tbl$allele1
    )
    if (!names_as_int) {
      site_names[[f]] <- scan$loci_tbl$marker.ID[new_loci]
    } else if (length(new_loci) > 0) {
      # the sidecar lines of the loci that created new sites, keyed by site
      scan_index <- sidecar_index(
        sidecars$scan, loci_sidecar_header, sidecar_block()
      )
      sidecar_gather(sidecars$scan, loci_sidecar_header,
        index_keys = scan_index$keys, index_offsets = scan_index$offsets,
        keys = new_loci, path_out = sidecars$sites,
        first_new_key = n_sites + 1, append = file.exists(sidecars$sites)
      )
    }
    unlink(sidecars$scan)
    n_sites <- n_sites + length(new_loci)
    rm(scan, new_loci)
  }
  sites <- merge_sites_result(merged_sites)
  if (sites$ref_mismatch$n > 0) {
    stop_ref_mismatch(sites, chromosome_names, files)
  }

  # the samples (in the order of the files)
  sample_id <- merge_vcfs_sample_names(sample_names, force_samples)
  n_indiv <- lengths(sample_names)

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
  used_names <- chromosome_names[chromosome_used]
  if (chromosomes_as_int) {
    # integer codes (see gen_tibble()), with their names in the sidecar
    chromosome_codes <- chromosome_int_codes_cpp(used_names)
    loci_chromosome <- chromosome_codes[
      chromosome_new_index[sites$chromosome[keep]]
    ]
    write_chromosome_sidecar_cpp(
      sidecars$chromosomes, chromosome_codes,
      used_names
    )
  } else {
    loci_chromosome <- chromosome_index_to_factor(
      chromosome_new_index[sites$chromosome[keep]],
      used_names
    )
  }
  loci_position <- sites$position[keep]
  if (names_as_int) {
    loci_name <- seq_along(keep)
    # the sidecar lines of the sites that are kept, in the order of the
    # merged loci (keyed by their integer names)
    sites_index <- sidecar_index(
      sidecars$sites, loci_sidecar_header, sidecar_block()
    )
    sidecar_gather(sidecars$sites, loci_sidecar_header,
      index_keys = sites_index$keys, index_offsets = sites_index$offsets,
      keys = keep, path_out = sidecars$loci,
      first_new_key = 1, append = FALSE
    )
    unlink(sidecars$sites)
  } else {
    # each locus takes the name it has in the first file where it is found
    # (missing IDs were filled with the chromosome and position when the
    # file was read)
    loci_name <- unlist(site_names)[keep]
    rm(site_names)
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
    missing_alleles = missing_alleles,
    chromosomes_as_int = chromosomes_as_int
  )

  # the column of the merged FBM (0-based) of each site (-1 if not kept)
  site_col <- rep(-1L, length(sites$status))
  site_col[keep] <- seq_along(keep) - 1L

  # second pass: the genotypes of each file, written straight into the
  # columns of its loci, in the rows of its samples (one file after the
  # other). With variants = "union", the samples of a file have missing
  # genotypes at the sites that are not in that file (or 0, with
  # missing_to_ref = TRUE, which is how an FBM starts), so we start with all
  # genotypes missing; with "intersection", all genotypes are written
  # (the columns of the loci of each file are only computed when its
  # genotypes are written)
  rm(sites)
  x_gt <- write_vcfs_gen_tibble(
    files = files,
    loci = loci,
    sample_id = sample_id,
    ploidy = unlist(file_ploidy),
    row_offset = c(0L, cumsum(n_indiv))[seq_len(n_files)],
    row_map = rep(list(integer(0)), n_files),
    file_cols = function(f) merge_sites_file_cols(merged_sites, f, site_col),
    file_ploidy = file_ploidy,
    fill_missing = variants == "union" && !missing_to_ref,
    backingfile = backingfile,
    n_cores = n_cores,
    quiet = quiet,
    phased = phased,
    file_samples = sample_names
  )
  sidecars_done <- TRUE
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
# files, listing the first few (the alleles of a site may come from several
# files, e.g. a missing ALT filled in by a later file, so they are described
# as those of the earlier files)
stop_ref_mismatch <- function(sites, chromosome_names, files) {
  mm <- sites$ref_mismatch
  alleles <- function(ref, alt) {
    paste0(ifelse(is.na(ref), ".", ref), "/", ifelse(is.na(alt), ".", alt))
  }
  examples <- vapply(seq_along(mm$site), function(i) {
    s <- mm$site[i]
    paste0(
      chromosome_names[sites$chromosome[s]], ":", sites$position[s], " (",
      alleles(mm$site_ref[i], mm$site_alt[i]), " in the earlier files, ",
      "first in ", basename(files[sites$first_file[s]]), "; ",
      alleles(mm$locus_ref[i], mm$locus_alt[i]), " in ",
      basename(files[mm$file[i]]), ")"
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
    "  loci with swapped alleles: ", fmt(sites$n_swapped), "\n",
    "  kept: ", fmt(length(keep)), " loci (", variants, ")"
  )
}

# The helpers below are shared by gt_merge_vcfs() and gt_concat_vcfs(), which
# read several VCF files into a single gen_tibble in two passes: a first pass
# reads the loci, samples and ploidy of each file (scan_vcf_files()), and a
# second pass writes the genotypes of each file straight into the backing file
# of the gen_tibble (write_vcfs_gen_tibble()).

# check the arguments shared by gt_merge_vcfs() and gt_concat_vcfs(),
# returning the normalised paths of the files
check_vcf_files_args <- function(files, valid_alleles, n_cores) {
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
  bad_n_cores <- !is.numeric(n_cores) || length(n_cores) != 1 ||
    !is.finite(n_cores) || n_cores < 1
  if (bad_n_cores) {
    stop("'n_cores' should be a positive integer")
  }
  files
}

# the backing file of the gen_tibble: by default, next to the first file, with
# a suffix (and a version number if the files already exist, as in
# gen_tibble())
default_vcfs_backingfile <- function(backingfile, files, suffix) {
  if (is.null(backingfile)) {
    backingfile <- paste0(
      sub("\\.(vcf\\.gz|vcf)$", "", files[1], ignore.case = TRUE),
      suffix
    )
  }
  change_duplicated_file_name(backingfile)
}

# the first pass: the loci of a file (no genotypes), with its samples and
# their ploidy (locus names are not read if we use integers: they are written
# to the sidecar loci_file instead, if it is not "")
scan_vcf_file <- function(file, names_as_int, n_cores, quiet,
                          loci_file = "") {
  if (!quiet) {
    message("reading the loci of ", file)
  }
  scan <- vcf_to_bk(file,
    bk_file = "", read_names = !names_as_int,
    n_threads = n_cores, read_genotypes = FALSE,
    loci_file = loci_file
  )
  if (length(scan$sample_names) == 0) {
    stop("the VCF ", file, " has no samples (or no variants)")
  }
  scan
}

# the first pass for all the files; with integer names, the loci of each
# file are added to the locus sidecar loci_file (keyed by their index in the
# concatenated files), if it is not ""
scan_vcf_files <- function(files, names_as_int, n_cores, quiet,
                           loci_file = "") {
  scans <- vector("list", length(files))
  n_loci <- 0
  for (f in seq_along(files)) {
    scan_file <- if (loci_file != "") paste0(loci_file, ".scan") else ""
    scans[[f]] <- scan_vcf_file(files[f], names_as_int, n_cores, quiet,
      loci_file = scan_file
    )
    if (loci_file != "") {
      # (all the lines of the file, keyed after those of the earlier files)
      sidecar_gather(scan_file, loci_sidecar_header,
        index_keys = integer(0), index_offsets = numeric(0),
        keys = integer(0), path_out = loci_file,
        first_new_key = n_loci + 1, append = f > 1
      )
      unlink(scan_file)
      n_loci <- n_loci + length(scans[[f]]$loci_tbl$physical.pos)
    }
  }
  scans
}

# the paths of the sidecars of a gen_tibble made from VCF files, next to its
# backing file (backingfile, without extension), and of their temporary
# files; sidecars left there by older files with the same name (which can not
# belong to the new data) are removed
vcfs_sidecar_paths <- function(backingfile) {
  bk_file <- paste0(backingfile, ".bk")
  remove_stale_sidecars(bk_file)
  loci <- sidecar_path_from_bk(bk_file, "loci")
  list(
    loci = loci,
    chromosomes = sidecar_path_from_bk(bk_file, "chromosomes"),
    # (the loci of a file being scanned, and the loci of the merged sites)
    scan = paste0(loci, ".scan"),
    sites = paste0(loci, ".sites")
  )
}

# the second pass: create the FBM of the gen_tibble, write the genotypes of
# each file into it, and create (and save) the gen_tibble. For each file,
# row_offset and row_map give the rows of its individuals (row_offset onwards,
# or the 0-based rows in row_map, if not empty), file_cols(f) gives a list
# with `col`, the 0-based column of each of its biallelic SNPs (-1 if not
# kept), and `swap`, whether the alleles of each SNP are swapped (or an empty
# vector), and file_ploidy the ploidy of its individuals. ploidy is the
# ploidy of the individuals of the gen_tibble. With phased = TRUE, phased
# genotypes are stored (see gen_tibble_phased()), and file_samples gives the
# names of the samples of each file (in the order of the file, for errors).
write_vcfs_gen_tibble <- function(files, loci, sample_id, ploidy, row_offset,
                                  row_map, file_cols, file_ploidy,
                                  fill_missing, backingfile, n_cores, quiet,
                                  phased = FALSE, file_samples = NULL) {
  if (anyNA(ploidy) || any(ploidy <= 0)) {
    stop("the ploidy of each individual should be a positive integer")
  }
  indiv_meta <- validate_indiv_meta(data.frame(id = sample_id))
  # the FBM: missing genotypes are its largest ploidy + 1 (all genotypes are
  # missing at the start if fill_missing, 0 otherwise, which is also 0|0 for
  # phased genotypes)
  max_ploidy <- max(ploidy)
  missing_value <- max_ploidy + 1
  if (phased) {
    # phased genotypes are only stored for diploid individuals (their missing
    # value is then 3, as in unphased data)
    stopifnot_phased_ploidy(ploidy, sample_id)
    code256 <- code_phased
  } else {
    code256 <- rep(NA_real_, 256)
    code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)
  }
  file_backed_matrix <- bigstatsr::FBM.code256(
    nrow = length(sample_id),
    ncol = nrow(loci),
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

  # the genotypes of each file, written straight into the columns of its
  # loci, in the rows of its samples
  for (f in seq_along(files)) {
    if (!quiet) {
      message("reading the genotypes of ", files[f])
    }
    this_file <- file_cols(f)
    vcf_merge_genotypes(files[f],
      BM = file_backed_matrix,
      row_offset = row_offset[f],
      row = row_map[[f]],
      col = this_file$col,
      swap = this_file$swap,
      ploidy = file_ploidy[[f]],
      missing_value = missing_value,
      n_threads = n_cores,
      phased_lut = if (phased) phased_genotype_byte else integer(0),
      sample_names = if (phased) file_samples[[f]] else character(0)
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
    # (dropping the gen_tibble, and the genotypes it was made from, first,
    # so that its backing file can be unmapped, see remove_fbm_files())
    gt_files <- gt_get_file_names(x_gt)
    rm(x_gt, indiv_meta)
    remove_fbm_files(gt_files)
    stop_genotypes_without_alleles()
  }
  gt_save_light(x_gt, quiet = quiet) # nolint
  x_gt
}
