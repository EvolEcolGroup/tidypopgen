# utility internal functions

# simple function to extract the extension of a file
# this avoids having to import tools just for tools::file_ext
file_ext <- function(x){
  utils::tail(unlist(strsplit(x,".",fixed = TRUE)),n=1)
}

# a simple tidier for dist matrices
# tidy.dist is deprecated, and often we have a full matrix rather than a dist object
tidy_dist_matrix<- function(mat){
  if (!inherits(mat,"matrix")){
    stop("mat should be a matrix")
  }
  mat_names<- colnames(mat)
  if (is.null(mat)){
    mat_names <- c(1:ncol(mat))
  }
  xy <- t(utils::combn(colnames(mat), 2))
  colnames(xy)<-c("item1","item2")
  xy %>% tibble::as_tibble() %>%
    dplyr::mutate(value=mat[xy])
}
