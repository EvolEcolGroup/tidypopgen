#' Cast a `gen_tibble` to a `genlight` object
#'
#' Convert an [gen_tibble] object to a [`adegenet::genlight`] object
#'
#' @param x the object to convert
#' @param ... additional arguments (currently unused)
#' @returns a `genlight` object
#' @rdname as_genlight
#' @export
#' @export
as_genlight <- function(x, ...) {
  new_gl <- new("genlight")
  new_gl@gen <- x$genotypes
  new_gl@n.loc <- nrow(show_loci(x))
  new_gl@ind.names <- x$id
  if ("population" %in% names(x)){
    new_gl@pop <- as.factor(x$population)
  }
  new_gl@loc.names<-show_loci_names(x)
  new_gl@loc.all <- paste0(show_loci(x)$allele_ref,"/",show_loci(x)$allele_alt) # TODO is this the right order?!?
  new_gl@chromosome <- as.factor(show_loci(x)$chromosome)
  new_gl@position <- as.integer(show_loci(x)$position)
  new_gl@ploidy <- unlist(lapply(x$genotypes,adegenet::ploidy))
  if (!validObject(new_gl)){
    stop("something went wrong when converting")
  } else{
    return(new_gl)
  }
}
