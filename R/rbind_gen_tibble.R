#' Combine two gen_tibbles
#'
#' This function combined two [gen_tibble]s. It subsets the loci and swaps ref
#' and alt alleles to make the two datasets compatible. The first object is
#' used as a
#' "reference" , and SNPs in the other dataset will be flipped and/or alleles
#' swapped as needed. SNPs that have different alleles in the two datasets (i.e.
#' triallelic) will also be dropped. There are also options (NOT default) to
#' attempt strand flipping to match alleles (often needed in human datasets from
#' different SNP chips), and remove ambiguous alleles (C/G and A/T) where the
#' correct strand can not be guessed.
#'
#' Loci are matched by name, unless `use_position = TRUE`. If the loci names of
#' either `gen_tibble` are integers (as created by `gen_tibble(names_as_int =
#' TRUE)`), loci are always matched by chromosome and position (see
#' `use_position`).
#'
#' rbind differs from merging data with plink, which swaps the order of allele1
#' and allele2 according to minor allele frequency when merging datasets. rbind
#' flips and/or swaps alleles according to the reference dataset, not according
#' to allele frequency.
#'
#' @param ... two [`gen_tibble`] objects. Note that this function can not take
#'   more objects, `rbind` has to be done sequentially for large sets of
#'   objects.
#' @param as_is deprecated, and ignored: `rbind` always subsets the loci and
#'   swaps alleles as needed (setting `as_is = TRUE` never stopped it, and now
#'   gives a warning).
#' @param use_position boolean of whether a combination of chromosome and
#'   position should be used for matching SNPs. By default, `rbind` uses the
#'   locus name, so this is set to FALSE. If the loci names of either
#'   `gen_tibble` are integers (as created by `gen_tibble(names_as_int =
#'   TRUE)`, which numbers the loci in the order in which they are found in the
#'   *VCF*), the names do not identify the same locus in the two datasets, so
#'   loci are always matched by position: `use_position` is then set to TRUE
#'   automatically (even if it was set to FALSE), with a message unless `quiet
#'   = TRUE`. In either case, the merged loci keep
#'   the names of the reference dataset. When using 'use_position=TRUE', make
#'   sure chromosomes are coded in the same way in both `gen_tibbles` (a mix of
#'   e.g. 'chr1', '1' or 'chromosome1' can be the reasons if an unexpectedly
#'   large number variants are dropped when merging).
#' @param flip_strand boolean on whether strand flipping should be checked to
#'   match the two datasets. If this is set to TRUE, ambiguous SNPs (i.e. A/T
#'   and C/G) will also be removed, and loci whose `allele_ref` is missing in
#'   both datasets will have their alleles swapped, so that the observed allele
#'   becomes `allele_ref`. It defaults to FALSE
#' @param quiet boolean whether to omit reporting to screen
#' @param backingfile the path and prefix of the files used to store the merged
#'   data (it will be a .RDS to store the `bigSNP` object and a .bk file as its
#'   backing file for the FBM)
#' @returns a [`gen_tibble`] with the merged data.
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Create a second gen_tibble to merge
#' test_indiv_meta <- data.frame(
#'   id = c("x", "y", "z"),
#'   population = c("pop1", "pop1", "pop2")
#' )
#' test_genotypes <- rbind(
#'   c(1, 1, 0, 1, 1, 0),
#'   c(2, 1, 0, 0, 0, 0),
#'   c(2, 2, 0, 0, 1, 1)
#' )
#' test_loci <- data.frame(
#'   name = paste0("rs", 1:6),
#'   chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
#'   position = as.integer(c(3, 5, 65, 343, 23, 456)),
#'   genetic_dist = as.double(rep(0, 6)),
#'   allele_ref = c("A", "T", "C", "G", "C", "T"),
#'   allele_alt = c("T", "C", NA, "C", "G", "A")
#' )
#'
#' test_gt <- gen_tibble(
#'   x = test_genotypes,
#'   loci = test_loci,
#'   indiv_meta = test_indiv_meta,
#'   valid_alleles = c("A", "T", "C", "G"),
#'   quiet = TRUE
#' )
#'
#' # Merge the datasets using rbind
#' merged_gt <- rbind(ref = example_gt, target = test_gt, flip_strand = TRUE)
#'
#' merged_gt
rbind.gen_tbl <- function(
  ...,
  as_is = FALSE,
  flip_strand = FALSE,
  use_position = FALSE,
  quiet = FALSE,
  backingfile = NULL
) {
  dots <- list(...)
  if (length(dots) != 2) {
    stop("rbind for gen_tibble can only take two tibbles at a time")
  }
  ref <- dots[[1]]
  target <- dots[[2]]
  # only bind diploid tibbles
  stopifnot_dip_pseudo(ref$genotypes)
  stopifnot_dip_pseudo(target$genotypes)
  any_pseudohaploid <- FALSE
  if (is_pseudohaploid(ref) || is_pseudohaploid(target)) {
    any_pseudohaploid <- TRUE
  }
  if (any(!is.na(match(ref$id, target$id)))) {
    stop(
      "The two gen_tibbles contain at least one individual with the",
      " same ID."
    )
  }
  # deprecation of as_is (which was never implemented: loci were always
  # harmonised)
  if (isTRUE(as_is)) {
    warning(
      "'as_is' is deprecated and ignored: rbind always harmonises the loci ",
      "(subsetting them and swapping alleles as needed)"
    )
  }
  # integer loci names can not be used to match loci, so match them by
  # position (rbind_dry_run() is then called with use_position = TRUE, so the
  # message is only given once)
  use_position <- rbind_use_position(
    ref_names = loci_names(ref),
    target_names = loci_names(target),
    use_position = use_position,
    quiet = quiet
  )
  # sort out paths
  if (is.null(backingfile)) {
    save_path <- dirname(attr(ref$genotypes, "fbm")$backingfile)
    backingfile <- tempfile("gt_merged_", tmpdir = save_path, fileext = "")
  }

  # check for duplicates within each dataset (if we use position, loci are
  # matched by a combination of chromosome and position, whose duplicates
  # are checked by rbind_dry_run(), so the names of the loci are left as
  # they are)
  if (!use_position) {
    if (any_duplicated_names(loci_names(target)) > 0) {
      stop(
        "The target gen_tibble contains duplicated loci. ",
        "Use 'find_duplicated_loci()' to identify and remove them."
      )
    }
    if (any_duplicated_names(loci_names(ref)) > 0) {
      stop(
        "The reference gen_tibble contains duplicated loci. ",
        "Use 'find_duplicated_loci()' to identify and remove them."
      )
    }
  }

  report <- rbind_dry_run(
    ref = ref,
    target = target,
    flip_strand = flip_strand,
    quiet = quiet,
    use_position = use_position
  )
  # now edit the gen_tibble objects
  ###########
  # we start with the ref object
  # we fix the missing alleles
  ## in the gt_table, in whichever allele column was missing
  for (allele_col in c("allele_alt", "allele_ref")) {
    id_missing <- which(report$ref$missing_allele_col %in% allele_col)
    attr(ref$genotypes, "loci")[[allele_col]][id_missing] <-
      report$ref$missing_allele[id_missing]
  }
  # swap the alleles of loci whose allele_ref is missing in both datasets, so
  # that the observed allele becomes allele_ref (only with flip_strand)
  ## in the gt_table (only if needed, to avoid copying the loci table)
  to_swap <- report$ref$to_swap
  if (any(to_swap)) {
    ref_loci <- attr(ref$genotypes, "loci")
    attr(ref$genotypes, "loci")$allele_ref[to_swap] <-
      ref_loci$allele_alt[to_swap]
    attr(ref$genotypes, "loci")$allele_alt[to_swap] <-
      ref_loci$allele_ref[to_swap]
  }
  # now create a new loci table (we'll use it later)
  ref_loci_order <- order(report$ref$new_id, na.last = NA)
  new_ref_loci_tbl <- show_loci(ref)[ref_loci_order, ]
  if (nrow(new_ref_loci_tbl) == 0) {
    stop("there are no loci in common between the two gen_tibbles")
  }
  ###########
  # now we move to the target object
  # we fix the missing alleles
  ## in the gt_table, in whichever allele column was missing
  for (allele_col in c("allele_alt", "allele_ref")) {
    id_missing <- which(report$target$missing_allele_col %in% allele_col)
    attr(target$genotypes, "loci")[[allele_col]][id_missing] <-
      report$target$missing_allele[id_missing]
  }
  # now flip the strands
  ## in the gt_table
  to_flip <- report$target$to_flip
  attr(target$genotypes, "loci")$allele_alt[to_flip] <-
    flip(attr(target$genotypes, "loci")$allele_alt[to_flip])
  attr(target$genotypes, "loci")$allele_ref[to_flip] <-
    flip(attr(target$genotypes, "loci")$allele_ref[to_flip])

  # now create a new loci table (we'll use it later)
  target_loci_order <- order(report$target$new_id, na.last = NA)
  new_target_loci_tbl <- show_loci(target)[target_loci_order, ]

  # now we create the merged FBM, with the individuals of ref followed by
  # those of target (it uses the code256 of ref)
  ref_rows <- vctrs::vec_data(ref$genotypes)
  target_rows <- vctrs::vec_data(target$genotypes)
  merged_fbm <- bigstatsr::FBM.code256(
    nrow = length(ref_rows) + length(target_rows),
    ncol = nrow(new_ref_loci_tbl),
    code = attr(ref$genotypes, "fbm")$code256,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )
  # and copy the raw genotypes of the loci in common directly into it (in
  # C++, without intermediate files), swapping the dosages of the loci whose
  # alleles are swapped (raw values 0 <-> 2, and imputed values 4 <-> 6): in
  # ref, loci whose allele_ref is missing in both datasets (only with
  # flip_strand), and in target, loci swapped to match ref
  swap_lut <- 0:255
  swap_lut[c(1, 3, 5, 7)] <- c(2, 0, 6, 4)
  gt_copy_fbm(
    BM_from = attr(ref$genotypes, "fbm"),
    BM_to = merged_fbm,
    rowInd = ref_rows,
    colInd = new_ref_loci_tbl$big_index,
    col_offset = 0,
    row_offset = 0,
    swap = report$ref$to_swap[ref_loci_order],
    swap_lut = swap_lut
  )
  gt_copy_fbm(
    BM_from = attr(target$genotypes, "fbm"),
    BM_to = merged_fbm,
    rowInd = target_rows,
    colInd = new_target_loci_tbl$big_index,
    col_offset = 0,
    row_offset = length(ref_rows),
    swap = report$target$to_swap[target_loci_order],
    swap_lut = swap_lut
  )
  # save the FBM (its .rds is used as fbm_file below)
  merged_fbm$save()

  # Now we need to create the gen_tibble
  # Make sure that the two tibbles have the same columns

  ref <- add_missing_cols(ref, target)
  names(ref$genotypes) <- ref$id
  target <- add_missing_cols(target, ref)
  names(target$genotypes) <- target$id
  merged_tbl <- rbind(
    ref %>% select(-dplyr::any_of("genotypes")),
    target %>% select(-dplyr::any_of("genotypes"))
  )
  # make sure that the genotypes vector points to the correct rows
  vctrs::vec_data(ref$genotypes)
  # and finally append the loci table
  indivs_with_big_names <- c(names(ref$genotypes), names(target$genotypes))
  # TODO check that this is the correct order!!!!
  new_ref_loci_tbl$big_index <- seq_len(nrow(new_ref_loci_tbl))
  # by default, this should just be a subset in the same order as the reference
  # TODO check that all individuals in tibble and bigsnp object are the same

  # the rds of the FBM (as for other gen_tibbles; its md5sum is quick to
  # compute, unlike that of the backing file)
  fbm_path <- merged_fbm$rds

  # fbm_ploidy covers the full FBM of each object, but only the rows in the
  # tibbles (in tibble order) were copied into the merged FBM
  indivs_ploidy <- c(
    attr(ref$genotypes, "fbm_ploidy")[vctrs::vec_data(ref$genotypes)],
    attr(target$genotypes, "fbm_ploidy")[vctrs::vec_data(target$genotypes)]
  )
  merged_tbl$genotypes <- new_vctrs_bigsnp(
    fbm_obj = merged_fbm,
    fbm_file = fbm_path,
    loci = new_ref_loci_tbl,
    indiv_id = indivs_with_big_names,
    ploidy = ifelse(any_pseudohaploid, -2L, 2L),
    fbm_ploidy = indivs_ploidy
  )

  merged_tibble <- tibble::new_tibble(
    merged_tbl,
    class = "gen_tbl"
  )
  gt_save(merged_tibble, quiet = quiet)
  return(merged_tibble)
}


add_missing_cols <- function(x, y) {
  missing_cols <- names(y)[!names(y) %in% names(x)]
  if (length(missing_cols) > 0) {
    for (col_i in missing_cols) {
      x[col_i] <- NA
    }
  }
  x
}
