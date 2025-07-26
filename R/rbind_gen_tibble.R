#' Combine two gen_tibbles
#'
#' This function combined two [gen_tibble]s. By defaults, it subsets the loci
#' and swaps ref and alt alleles to make the two datasets compatible (this
#' behaviour can be switched off with `as_is`). The first object is used as a
#' "reference" , and SNPs in the other dataset will be flipped and/or alleles
#' swapped as needed. SNPs that have different alleles in the two datasets (i.e.
#' triallelic) will also be dropped. There are also options (NOT default) to
#' attempt strand flipping to match alleles (often needed in human datasets from
#' different SNP chips), and remove ambiguous alleles (C/G and A/T) where the
#' correct strand can not be guessed.
#'
#' rbind differs from merging data with plink, which swaps the order of allele1
#' and allele2 according to minor allele frequency when merging datasets. rbind
#' flips and/or swaps alleles according to the reference dataset, not according
#' to allele frequency.
#'
#' @param ... two [`gen_tibble`] objects. Note that this function can not take
#'   more objects, `rbind` has to be done sequentially for large sets of
#'   objects.
#' @param as_is boolean determining whether the loci should be left as they are
#'   before merging. If FALSE (the defaults), `rbind` will attempt to subset and
#'   swap alleles as needed.
#' @param use_position boolean of whether a combination of chromosome and
#'   position should be used for matching SNPs. By default, `rbind` uses the
#'   locus name, so this is set to FALSE. When using 'use_position=TRUE', make
#'   sure chromosomes are coded in the same way in both `gen_tibbles` (a mix of
#'   e.g. 'chr1', '1' or 'chromosome1' can be the reasons if an unexpectedly
#'   large number variants are dropped when merging).
#' @param flip_strand boolean on whether strand flipping should be checked to
#'   match the two datasets. If this is set to TRUE, ambiguous SNPs (i.e. A/T
#'   and C/G) will also be removed. It defaults to FALSE
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
    backingfile = NULL) {
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
  if (any(!is.na(match(
    attr(ref$genotypes, "bigsnp")$fam$sample.ID,
    attr(target$genotypes, "bigsnp")$fam$sample.ID
  )))) {
    stop(
      "The two bigsnp objects contain at least one individual with the",
      " same ID."
    )
  }
  if (!quiet) {
    if (as_is) {
      if (flip_strand) {
        stop(paste(
          "if 'as_is' is set to TRUE, 'flip_strand' and",
          "'remove_ambiguous' have to be FALSE"
        ))
      }
      message(
        "rbind will not attempt to harmonise the loci (e.g. flip strand,\n",
        " reorder or subset) set 'as_is' to FALSE to subset to loci in common\n"
      )
    }
  }
  # sort out paths
  if (is.null(backingfile)) {
    save_path <- dirname(attr(ref$genotypes, "bigsnp")$genotypes$backingfile)
    backingfile <- tempfile("gt_merged_", tmpdir = save_path, fileext = "")
  }

  # if we use position, we update the names of the loci
  if (use_position) {
    show_loci(ref) <-
      show_loci(ref) %>%
      mutate(
        name_old = .data$name,
        name = paste(.data$chromosome, .data$position, sep = "_")
      )
    show_loci(target) <-
      show_loci(target) %>%
      mutate(
        name_old = .data$name,
        name = paste(.data$chromosome, .data$position, sep = "_")
      )
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
  ## in the gt_table
  id_missing <- which(!is.na(report$ref$missing_allele))
  attr(ref$genotypes, "loci")$allele_alt[id_missing] <-
    report$ref$missing_allele[id_missing]
  # and in  the bigSNP object
  attr(ref$genotypes, "bigsnp")$map$allele1[attr(
    ref$genotypes,
    "loci"
  )$big_index[id_missing]] <- # nolint
    report$ref$missing_allele[id_missing]
  # now create a new loci table (we'll use it later)
  new_ref_loci_tbl <- show_loci(ref)[order(report$ref$new_id, na.last = NA), ]
  # now we subset the SNP object
  ## in the snp object
  if (nrow(new_ref_loci_tbl) == 0) {
    stop("there are no loci in common between the two gen_tibbles")
  }
  ref_snp <- subset_bigSNP(
    attr(ref$genotypes, "bigsnp"),
    loci_indices = new_ref_loci_tbl$big_index,
    indiv_indices = vctrs::vec_data(ref$genotypes)
  )
  ###########
  # now we move to the target object
  # we fix the missing alleles
  ## in the gt_table
  id_missing <- which(!is.na(report$target$missing_allele))
  attr(target$genotypes, "loci")$allele_alt[id_missing] <-
    report$target$missing_allele[id_missing]
  # and in  the bigSNP object
  attr(target$genotypes, "bigsnp")$map$allele1[attr(
    target$genotypes,
    "loci"
  )$big_index[id_missing]] <- # nolint
    report$target$missing_allele[id_missing]
  # now flip the strands
  ## in the gt_table

  to_flip <- report$target$to_flip
  attr(target$genotypes, "loci")$allele_alt[to_flip] <-
    flip(attr(target$genotypes, "loci")$allele_alt[to_flip])
  attr(target$genotypes, "loci")$allele_ref[to_flip] <-
    flip(attr(target$genotypes, "loci")$allele_ref[to_flip])
  to_flip_big_index <- show_loci(target)$big_index[to_flip]
  attr(target$genotypes, "bigsnp")$map$allele1[to_flip_big_index] <-
    flip(attr(target$genotypes, "bigsnp")$map$allele1[to_flip_big_index])
  attr(target$genotypes, "bigsnp")$map$allele2[to_flip_big_index] <-
    flip(attr(target$genotypes, "bigsnp")$map$allele2[to_flip_big_index])

  # now create a new loci table (we'll use it later)
  new_target_loci_tbl <-
    show_loci(target)[order(report$target$new_id, na.last = NA), ]
  # now we subset the SNP object
  ## in the snp object
  target_snp <- subset_bigSNP(
    attr(target$genotypes, "bigsnp"),
    loci_indices = new_target_loci_tbl$big_index,
    indiv_indices = vctrs::vec_data(target$genotypes),
    swap_indices = show_loci(target)$big_index[report$target$to_swap]
  )

  # this check does not work if we change the names of markers
  #  if (!identical(target_snp$map$marker.ID,ref_snp$map$marker.ID) |
  #      !identical(target_snp$map$allele2,ref_snp$map$allele2)) {
  #    stop("something went wrong when subsetting and
  #          reordering the two bigSNP objects")
  #  }

  # now we need to merge the two FBMs
  # we start by transposing them, so that they just need to be appended
  t_ref_fbm <- bigstatsr::big_transpose(ref_snp$genotypes)
  t_target_fbm <- bigstatsr::big_transpose(target_snp$genotypes)
  # append the two files
  append_success <- file.append(t_ref_fbm$backingfile, t_target_fbm$backingfile) # nolint
  # and amend the new number of columns
  t_ref_fbm$ncol <- t_ref_fbm$ncol + t_target_fbm$ncol
  # now flip the file around
  merged_fbm <- bigstatsr::big_transpose(t_ref_fbm, backingfile = backingfile)
  # TODO this should be written in the directory of interest

  # Make sure that the two fam tibbles have the same columns
  ref_snp$fam <- add_missing_cols(ref_snp$fam, target_snp$fam)
  target_snp$fam <- add_missing_cols(target_snp$fam, ref_snp$fam)

  # now create a bigsnp object
  merged_snp <- structure(
    list(
      genotypes = merged_fbm,
      fam = rbind(ref_snp$fam, target_snp$fam),
      map = ref_snp$map
    ),
    class = "bigSNP"
  )
  merged_rds <- paste0(backingfile, ".rds")
  saveRDS(merged_snp, merged_rds)
  # Now we need to create the gen_tibble
  # Make sure that the two tibbles have the same columns
  ref <- add_missing_cols(ref, target)
  target <- add_missing_cols(target, ref)
  merged_tbl <- rbind(
    ref %>% select(-dplyr::any_of("genotypes")),
    target %>% select(-dplyr::any_of("genotypes"))
  )
  # make sure that the genotypes vector points to the correct rows
  vctrs::vec_data(ref$genotypes)
  # and finally append the loci table
  indivs_with_big_names <- c(names(ref$genotypes), names(target$genotypes))
  # new_ref_loci_tbl$big_index<-
  #               match(new_ref_loci_tbl$name,merged_snp$map$marker.ID) #nolint
  # TODO check that this is the correct order!!!!
  new_ref_loci_tbl$big_index <- seq_len(nrow(new_ref_loci_tbl))
  # by default, this should just be a subset in the same order as the reference
  # TODO check that all individuals in tibble and bigsnp object are the same
  merged_tbl$genotypes <- vctrs::new_vctr(
    match(
      indivs_with_big_names,
      merged_snp$fam$sample.ID
    ),
    # TODO check that this is the correct order!!!!
    bigsnp = merged_snp,
    bigsnp_file = merged_rds,
    bigsnp_md5sum = tools::md5sum(merged_rds),
    loci = new_ref_loci_tbl,
    names = indivs_with_big_names,
    # TODO currently set to only work for diploids or pseudohaploids
    ploidy = ifelse(any_pseudohaploid, -2L, 2L),
    class = "vctrs_bigSNP"
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
