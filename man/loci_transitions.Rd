% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loci_transitions.R
\name{loci_transitions}
\alias{loci_transitions}
\alias{loci_transitions.tbl_df}
\alias{loci_transitions.vctrs_bigSNP}
\alias{loci_transitions.grouped_df}
\title{Find transitions}
\usage{
loci_transitions(.x, ...)

\method{loci_transitions}{tbl_df}(.x, ...)

\method{loci_transitions}{vctrs_bigSNP}(.x, ...)

\method{loci_transitions}{grouped_df}(.x, ...)
}
\arguments{
\item{.x}{a vector of class \code{vctrs_bigSNP} (usually the \code{genotype} column of
a \code{\link{gen_tibble}} object),
or a \code{\link{gen_tibble}}.}

\item{...}{other arguments passed to specific methods.}
}
\value{
a logical vector defining which loci are transitions
}
\description{
Use the loci table to define which loci are transitions
}
