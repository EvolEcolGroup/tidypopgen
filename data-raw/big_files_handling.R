library(bigstatsr)
X <- FBM.code256(10, 10, code=bigsnpr::CODE_012)
line <- sample(c(0,1,2,4),10,replace=TRUE)
X[1,] <-line
X$is_saved
X$save()

X$save()

first_file <- X$backingfile

rm(X)
X3 <- bigstatsr::big_attach(first_file)


X2 <- bigstatsr::big_copy(X)
new_file <- X2$backingfile

rm(X,X2)


