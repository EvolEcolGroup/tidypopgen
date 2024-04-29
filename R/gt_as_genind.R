#' Convert a `gen_tibble` to a `genind` object from `adegenet`
#'
#' This function converts a `gen_tibble` to a `genind` object from `adegenet`
#'
#' @param x a [`gen_tibble`], with population coded as 'population'
#' @returns a `genind` object
#' @export
#'

gt_as_genind <- function(x){
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop(
      "to use this function, first install package 'adegenet' with\n",
      "install.packages('adegenet')")
  }
  df_for_genind <- show_genotypes(x)
  df_for_genind [df_for_genind ==0]<-"11"
  df_for_genind [df_for_genind ==1]<-"12"
  df_for_genind [df_for_genind ==2]<-"22"
  test_genind <- adegenet::df2genind(X = df_for_genind,
                                     ind.names = x$id,
                                     pop = x$population,
                                     ncode=1,
                                     loc.names = show_loci_names(x))
}

