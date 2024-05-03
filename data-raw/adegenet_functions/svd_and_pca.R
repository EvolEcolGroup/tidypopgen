# https://genomicsclass.github.io/book/pages/pca_svd.html
# The elements of D are formed by taking the sum of the squares of the principal
# components but not dividing by the sample size.

# We need to center the matrix for the svd to make sense, and ideally scale
# if we scale, then the


# compare svd to pca
# create a random dataset (10 ind and 20 features)
a <- FBM(10, 20)
big_apply(a, a.FUN = function(X, ind) {
  X[, ind] <- rnorm(nrow(X) * length(ind))
  NULL
}, a.combine = 'c')

a_mat <- a[]


a_prcomp <- prcomp(a_mat,center=TRUE, scale=FALSE)
a_cent <- sweep(a_mat, 2, colMeans(a_mat), "-")
a_svd <- svd(a_cent)

# U*D in svd is x in prcomp
all.equal(
  a_prcomp$x,
  sweep(a_svd$u, 2, a_svd$d, '*'),
  check.attributes = FALSE)

# V in svd is rotation in prcomp
all.equal(
  a_prcomp$rotation,
  a_svd$v,
  check.attributes = FALSE)

# eigen values
## eigen
all.equal(
  a_prcomp$sdev^2,
  a_svd$d^2/(nrow(a) - 1),
  check.attributes = FALSE)

# for an centered matrix (but not scaled), the sum of the eigen vectors is equal
#to the sum of variances
all.equal(sum(a_prcomp$sdev^2), sum(apply(a_mat,2,var)))




