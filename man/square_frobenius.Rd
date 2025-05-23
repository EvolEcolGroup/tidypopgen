% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/square_frobenius.R
\name{square_frobenius}
\alias{square_frobenius}
\title{Compute the square of the Frobenius norm of a matrix}
\usage{
square_frobenius(
  X,
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  center = rep(0, length(ind.col)),
  scale = rep(1, length(ind.col))
)
}
\arguments{
\item{X}{An FBM matrix}

\item{ind.row}{A vector of row indices}

\item{ind.col}{A vector of column indices}

\item{center}{A vector of centering values (i.e the means of the genotype
counts)}

\item{scale}{A vector of scaling values (i.e the standard deviations of the
genotype counts)}
}
\value{
The square of the Frobenius norm of the matrix
}
\description{
This function computes the square of the Frobenius norm of a matrix, which is
the sum of the squares of the matrix elements, which provides a measure of
of the total variance of the matrix.
The code used here was outlined in an issue of bigstatsr by @privefl:
https://github.com/privefl/bigstatsr/issues/83
}
\keyword{internal}
