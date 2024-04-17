#' Convert vcf to bed
#'
#' This function quickly converts a smallish vcf to bed. The vcf needs to be
#' read into memory and manipulated in memory, so there are some serious limits
#' to this function. PLINK is a better option for converting vcf to BED, but
#' for small datasets, it should do the job.
#' @param vcf_path the path to the vcf
#' @param bed_path the path where the bed file should be saved. If left to null,
#' the vcf path will be used by switching the file type from .vcf to .bed
#' @param quiet if TRUE, don't output information to screen
#' @param ... parameters to be passed to [vcfR::read.vcfR]
#' @returns the file path where the BED file was saved.
#' @export

gt_vcf_to_bed <- read_vcf <- function(vcf_path, bed_path=NULL, quiet = TRUE, ...) {

  if (requireNamespace("vcfR", quietly = TRUE)) {
    .x <- vcfR::read.vcfR(file = vcf_path, verbose = !quiet, ...)
    # subset to biallelic loci only
    bi <- vcfR::is.biallelic(.x)
    if(sum(!bi) > 0){
      msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
      msg <- c(msg, "\n", paste("Objects of class gen_tibble only support loci with two alleles."))
      msg <- c(msg, "\n", paste(sum(!bi), 'loci will be omitted from the gen_tibble object.'))
      warning(msg)
      .x <- .x[bi,]
    }

    # fill in any missing IDs of loci
    .x <- vcfR::addID(.x)

    # create loci table
    loci <- tibble(name = vcfR::getID(.x),
                   chromosome = vcfR::getCHROM(.x),
                   position = vcfR::getPOS(.x),
                   allele_ref = vcfR::getREF(.x),
                   allele_alt = vcfR::getALT(.x))

    .x <- vcfR::extract.gt(.x)
    .x[.x=="0|0"] <- 0
    .x[.x=="0|1"] <- 1
    .x[.x=="1|0"] <- 1
    .x[.x=="1|1"] <- 2
    .x[.x=="0/0"] <- 0
    .x[.x=="0/1"] <- 1
    .x[.x=="1/0"] <- 1
    .x[.x=="1/1"] <- 2
    # remove warnings about NAs
    suppressWarnings(storage.mode(.x) <- "integer")

    indiv_meta <- tibble(id = colnames(.x), population = NA)

    if (is.null(bed_path)){
      bed_path <- sub_vcf(vcf_path)
    } else if (file_ext(bed_path)=="bed"){
      bed_path <- bigsnpr::sub_bed(bed_path,"")
    }
    bed_path <- gt_write_bed_from_dfs(indiv_meta= indiv_meta,
                          genotypes = t(.x),
                          loci = loci,
                          path_out = bed_path)

    return(bed_path)
  } else {
    stop(
      "to convert from vcfR objects, first install package 'vcfR' with\n",
      "install.packages('vcfR')"
    )
  }
}

# this adapts bigsnpr::sub_bed to vcf
sub_vcf <- function (path, replacement = "", stop_if_not_ext = TRUE)
{
  pattern <- "\\.vcf$"
  if (!grepl(pattern, path))
    stop("Path '%s' must have 'vcf' extension.", path)
  if (stop_if_not_ext && (nchar(replacement) > 0) && (substr(replacement,
                                                             1, 1) != "."))
    stop("Replacement must be an extension starting with '.' if provided.")
  sub(pattern, replacement, path)
}
