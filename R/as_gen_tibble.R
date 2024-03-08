#' Cast to a `gen_tibble`
#'
#' Convert an object to a `gen_tibble`
#'
#' @param .x the object to convert
#' @param ... additional arguments to pass to a given method
#' @returns a `gen_tbl` object
#' @rdname as_gen_tibble
#' @export


#' @export
as_gen_tibble <- function(.x, ...) {
  UseMethod("as_gen_tibble", .x)
}

#################################################################
## method for genlight
#################################################################


#' @rdname as_gen_tibble
#' @param ignore_null_slots boolean on whether null slots should be ignored
#' (and thus filled with default values), or whether an error will be thrown
#' @export

as_gen_tibble.genlight <- function(.x, ..., ignore_null_slots = TRUE){
  n_ind <- adegenet::nInd(.x)
  n_loci <- adegenet::nLoc(.x)
  ind_meta_list <- list(id = fix_null_slot(.x,"ind.names",
                                           paste0("ind",seq_len(n_ind)),
                                           ignore_null_slots),
                        population = fix_null_slot(.x,"pop",
                                                   rep("pop1",n_ind),
                                                   ignore_null_slots))
  loci <- tibble(name = fix_null_slot(.x,"loc.names",
                                      paste0("loc",seq_len(n_loci)),
                                      ignore_null_slots),
                 chromosome = fix_null_slot(.x,"chromosome",
                                            rep("1",n_loci),
                                            ignore_null_slots),
                 position = fix_null_slot(.x,"position",
                                          seq_len(n_loci),
                                          ignore_null_slots))
  # now deal with alleles
  if (!is.null(.x@loc.all)){
    loci$allele_ref <- substr(.x@loc.all,1,1)
    loci$allele_alt <- substr(.x@loc.all,3,3)
  } else {
    loci$allele_ref <- fix_null_slot(.x,"loc.all",
                                    rep("x",n_loci),
                                    ignore_null_slots)
    loci$allele_alt <- fix_null_slot(.x,"loc.all",
                                    rep("y",n_loci),
                                    ignore_null_slots)
  }
  ind_meta_list$genotypes <- .x$gen
  attr(ind_meta_list$genotypes,"loci")<-loci

  tibble::new_tibble(
    ind_meta_list,
    class = "gen_tbl"
  )
}

# quick function to replace null slots with default replacement value or throw error
fix_null_slot <- function(object, slot_string, replacement, ignore_null_slots=TRUE){
  if (is.null(methods::slot(object,slot_string))){
    if (ignore_null_slots){
      return(replacement)
    } else {
      stop(slot_string, " is null and `ignore_null_slots`=FALSE")
    }
  } else {
    return(methods::slot(object,slot_string))
  }
}

#################################################################
## method for vcfR
#################################################################

#' @rdname as_gen_tibble
#' @param n_cores currently unused, for future parallelisation
#' @export

as_gen_tibble.vcfR<- function(.x, n_cores=1,...){
  if (requireNamespace("vcfR", quietly = TRUE)) {
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

    ind_meta <- tibble(id = colnames(.x), population = NA)

    .x <- gen_tibble(ind_meta= ind_meta, genotypes = .x, loci = loci)

    return(.x)
  } else {
    stop(
      "to convert from vcfR objects, first install package 'vcfR' with\n",
      "install.packages('vcfR')"
    )
  }
}
