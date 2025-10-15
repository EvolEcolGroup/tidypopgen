#' Convert vcf to FBM using vcfR as a parser.
#'
#' Convert a vcf file to a Filebacked Big Matrix (FBM) object. This should work
#' even for large vcf files that would not fit in memory.
#'
#' @param vcf_path the path to the vcf
#' @param chunk_size the chunk size to use on the vcf when loading the file
#' @param backingfile the name of the file to use as the backing file
#' @param valid_alleles a character vector of valid alleles. Default is c("A",
#'   "T", "C", "G").
#' @param missing_alleles a character vector of alleles to be treated as
#'   missing. Default is c("0", ".").
#' @param allow_duplicates whether to allow duplicated loci (same chromosome and
#'   position) or duplicated locus names. Default is FALSE.
#' @param quiet whether to print messages
#' @param ... further arguments to be passed to [vcfR::read.vcfR()]. Do not pass
#'   nrows, skip, verbose, or convertNA; these are controlled internally.
#' @returns an object of the class `gen_tbl`.
#' @keywords internal
#' @noRd
# nolint start
vcf_to_fbm_vcfR <- function(
    # nolint end
    vcf_path,
    chunk_size = NULL,
    backingfile = NULL,
    valid_alleles = c("A", "T", "C", "G"),
    missing_alleles = c("0", "."),
    allow_duplicates = FALSE,
    quiet = FALSE,
    ...) {
  dots <- list(...)
  forbidden <- intersect(
    names(dots),
    c("nrows", "skip", "verbose", "convertNA")
  )
  if (length(forbidden)) {
    stop("Unsupported via ...: ", paste(forbidden, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(backingfile)) {
    backingfile <- vcf_path
    backingfile <- sub("\\.vcf.gz$", "", backingfile)
    backingfile <- sub("\\.vcf$", "", backingfile)
  }
  if (file_ext(backingfile) == "bk") {
    backingfile <- bigstatsr::sub_bk(backingfile, "")
  }

  # count the variants in the file
  no_variants <- count_vcf_variants(vcf_path)
  # count individuals in the file
  no_individuals <- count_vcf_individuals(vcf_path)
  # if chunk is null, get the best guess of an efficient approach
  if (is.null(chunk_size)) {
    chunk_size <- bigstatsr::block_size(no_individuals)
  }

  # set up chunks
  chunks_vec <- c(
    rep(chunk_size, floor(no_variants / chunk_size)),
    no_variants %% chunk_size
  )
  # remove any 0 length chunks (there could be one at the end)
  chunks_vec <- chunks_vec[chunks_vec > 0L]

  # figure out ploidy from the first marker
  temp_vcf <- vcfR::read.vcfR(
    vcf_path,
    nrows = 1,
    verbose = !quiet,
    convertNA = FALSE,
    ...
  )
  temp_gt <- vcfR::extract.gt(temp_vcf, convertNA = FALSE)
  ploidy <- unname(apply(temp_gt, 2, get_ploidy))


  # check that there are no missing values in ploidy vector
  if (any(is.na(ploidy))) {
    stop("'ploidy' can not contain NAs")
  }
  # check all values > 0
  if (any(ploidy <= 0)) {
    stop(
      "the vector of individual ploidies ('ploidy') must contain ",
      "positive integers"
    )
  }
  fbm_ploidy <- ploidy
  max_ploidy <- max(ploidy)


  # set up codes for the appropriate ploidy level
  code256 <- rep(NA_real_, 256)
  code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)

  indiv_meta <- list(
    id = colnames(temp_gt)
  )

  loci <- tibble::tibble(
    name = NULL,
    chromosome = NULL,
    position = NULL,
    genetic_dist = NULL,
    allele_ref = NULL,
    allele_alt = NULL
  )

  # create the file backed matrix
  file_backed_matrix <- bigstatsr::FBM.code256(
    nrow = no_individuals,
    ncol = 0,
    code = code256,
    backingfile = backingfile
  )

  for (i in seq_along(chunks_vec)) {
    temp_vcf <- vcfR::read.vcfR(
      vcf_path,
      nrows = chunks_vec[i],
      skip = sum(chunks_vec[seq_len(i - 1)]),
      verbose = !quiet,
      ...
    )
    # filter any marker that is not biallelic
    bi <- vcfR::is.biallelic(temp_vcf)
    gt <- vcfR::extract.gt(temp_vcf, convertNA = FALSE)
    gt <- gt[bi, , drop = FALSE]
    if (nrow(gt) > 1) {
      # @TODO we could parallelise here
      gt <- t(apply(gt, 2, poly_indiv_dosage, max_ploidy = max_ploidy))
    } else if (nrow(gt) == 1) {
      # if we only have one marker
      gt <-
        matrix(
          apply(gt, 2, poly_indiv_dosage, max_ploidy = max_ploidy),
          ncol = 1
        )
    } else {
      next
    }
    # expand the file backed matrix according to the size of the gt matrix
    # get current size
    index_start <- dim(file_backed_matrix)[2] + 1
    # add the new columns
    file_backed_matrix$add_columns(ncol(gt))
    # fill them in
    file_backed_matrix[
      ,
      index_start:(index_start + ncol(gt) - 1)
    ] <- gt

    # add metadata
    # empty metadata
    temp_vcf <- vcfR::addID(temp_vcf)

    # create loci table
    loci <- rbind(
      loci,
      tibble(
        # remove names as it does have ID as a name
        name = unname(vcfR::getID(temp_vcf)[bi]),
        chromosome = unname(vcfR::getCHROM(temp_vcf)[bi]),
        position = vcfR::getPOS(temp_vcf)[bi],
        genetic_dist = 0,
        allele_ref = unname(vcfR::getREF(temp_vcf)[bi]),
        allele_alt = unname(vcfR::getALT(temp_vcf)[bi])
      )
    )
  }
  # validate the loci
  loci <- validate_loci(loci,
    check_alphabet = TRUE,
    harmonise_loci = TRUE,
    check_duplicates = TRUE,
    allow_duplicates = allow_duplicates,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles
  )
  # validate individuals
  indiv_meta <- validate_indiv_meta(as.data.frame(indiv_meta))
  # construct path
  fbm_path <- bigstatsr::sub_bk(file_backed_matrix$backingfile, ".rds")

  # create genotypes column
  indiv_meta$genotypes <- new_vctrs_bigsnp(
    fbm_obj = file_backed_matrix,
    fbm_file = fbm_path,
    loci = loci,
    indiv_id = indiv_meta$id,
    ploidy = max_ploidy,
    fbm_ploidy = fbm_ploidy
  )

  new_gen_tbl <- tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )
  return(new_gen_tbl)
}

# get ploidy for a given individual
get_ploidy <- function(x) {
  max(sapply(strsplit(x, "[/|]"), function(x) length(x)))
}

# get dosages for all the genotypes of an individual
# (x is a vector of genotypes in the standard vcf format)
poly_indiv_dosage <- function(x, max_ploidy) {
  sapply(strsplit(x, "[/|]"), poly_genotype_dosage, max_ploidy)
}

# get dosages for a genotype x and return them as raw for inclusion in the FBM
poly_genotype_dosage <- function(x, max_ploidy) {
  if (!is.na(x[1]) && x[1] != ".") {
    dosage <- sum(as.numeric(x))
    if (dosage < (max_ploidy + 1)) {
      return(as.raw(dosage)) # nolint
    } else {
      stop(paste(
        "a genotype has more than max_ploidy alleles. We estimate",
        "max_ploidy from the first variant in the vcf file, make",
        "sure that variant is representative of ploidy (e.g. it is",
        "not on a sex chromosome)."
      ))
    }
  } else {
    return(as.raw(max_ploidy + 1)) # nolint
  }
}
