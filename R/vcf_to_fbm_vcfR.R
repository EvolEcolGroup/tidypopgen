#' Convert vcf to FBM using vcfR as a parser.
#'
#' Convert a vcf file to a Filebacked Big Matrix (FBM) object.
#' This should work even for large vcf files that would not fit in memory.
#'
#' @param vcf_path the path to the vcf
#' @param chunk_size the chunk size to use on the vcf when loading the file
#' @param backingfile the name of the file to use as the backing file
#' @return path to the resulting rds file as class bigSNP.
#' @keywords internal
#' @noRd
# nolint start
vcf_to_fbm_vcfR <- function(
    # nolint end
    vcf_path,
    chunk_size = NULL,
    backingfile = NULL,
    quiet = FALSE) {
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

  # figure out ploidy from the first marker
  temp_vcf <- vcfR::read.vcfR(
    vcf_path,
    nrow = 1,
    verbose = !quiet,
    convertNA = FALSE
  )
  temp_gt <- vcfR::extract.gt(temp_vcf, convertNA = FALSE)
  ploidy <- apply(temp_gt, 2, get_ploidy)
  if (any(is.na(ploidy))) {
    stop("error whilst determining ploidy")
  }
  max_ploidy <- max(ploidy)

  # set up codes for the appropriate ploidy level
  code256 <- rep(NA_real_, 256)
  code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)

  # metadata
  fam <- tibble(
    family.ID = colnames(temp_gt),
    sample.ID = colnames(temp_gt),
    paternal.ID = 0,
    maternal.ID = 0,
    sex = 0,
    affection = -9,
    ploidy = unname(ploidy)
  )

  loci <- tibble(
    chromosome = NULL,
    marker.id = NULL,
    genetic.dist = NULL,
    physical.pos = NULL,
    allele1 = NULL,
    allele2 = NULL
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
      nrow = chunks_vec[i],
      skip = sum(chunks_vec[seq_len(i - 1)]),
      verbose = !quiet
    )
    # filter any marker that is not biallelic
    bi <- vcfR::is.biallelic(temp_vcf)
    gt <- vcfR::extract.gt(temp_vcf)
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
        chromosome = unname(vcfR::getCHROM(temp_vcf)[bi]),
        # remove names as it does have ID as a name
        marker.ID = unname(vcfR::getID(temp_vcf)[bi]),
        genetic.dist = 0,
        physical.pos = vcfR::getPOS(temp_vcf)[bi],
        allele1 = unname(vcfR::getALT(temp_vcf)[bi]),
        allele2 = unname(vcfR::getREF(temp_vcf)[bi])
      )
    )
  }
  # save it
  file_backed_matrix$save()

  bigsnp_obj <- structure(
    list(
      genotypes = file_backed_matrix,
      fam = fam,
      map = loci
    ),
    class = "bigSNP"
  )

  bigsnp_obj <- bigsnpr::snp_save(bigsnp_obj)
  # and return the path to the rds
  bigsnp_obj$genotypes$rds
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
