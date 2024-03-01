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
