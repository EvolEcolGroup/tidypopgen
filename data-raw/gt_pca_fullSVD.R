#' PCA for `gen_tibble` objects by full SVD
#'
#' This function performs Principal Component Analysis on a `gen_tibble`,
#' by computing the full SVD through eigen decomposition of the covariance matrix.
#'  This function is based on code posted to
#' StackOverflow[https://stackoverflow.com/questions/46253537/computing-the-null-space-of-a-bigmatrix-in-r/46380540#46380540]
#'  by @privefl
#'  @export


big_fullSVD <- function (X,
                         fun.scaling = big_scale(center = FALSE, scale = FALSE),
                         ind.row = rows_along(X),
                         ind.col = cols_along(X),
                         k = 10){

  # Compute t(a) * a
  K <- big_crossprodSelf(X, fun.scaling)

  # Get v and d where a = u * d * t(v) the SVD of a
  eig <- eigen(K[])
  v <- eig$vectors
  d <- sqrt(eig$values)

  # we should esimate the rank here and crop the vectors and values as in the pca for adegenet


  means <- attr(K, "center")
  sds   <- attr(K, "scale")

  # Get u if you need it. It will be of the same size of u
  # so that I store it as a FBM.
  u <- FBM(nrow(X), ncol(X))
  big_apply(u, a.FUN = function(X, ind, a, v, d) {
    X[ind, ] <- sweep(a[ind, ] %*% v, 2, d, "/")
    NULL
  }, a.combine = 'c', block.size = 50e3, ind = rows_along(u),
  a = X, v = v, d = d)

  res <-list(d=d, u=u, v=v, center = means, scale= sds)
  structure(res, class = "big_SVD")
}

