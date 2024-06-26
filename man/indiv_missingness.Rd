% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indiv_missingness.R
\name{indiv_missingness}
\alias{indiv_missingness}
\alias{indiv_missingness.tbl_df}
\alias{indiv_missingness.vctrs_bigSNP}
\alias{indiv_missingness.grouped_df}
\title{Estimate individual missingness}
\usage{
indiv_missingness(.x, as_counts = FALSE, ...)

\method{indiv_missingness}{tbl_df}(.x, as_counts = FALSE, ...)

\method{indiv_missingness}{vctrs_bigSNP}(.x, as_counts = FALSE, ...)

\method{indiv_missingness}{grouped_df}(.x, as_counts = FALSE, ...)
}
\arguments{
\item{.x}{a vector of class \code{vctrs_bigSNP} (usually the \code{genotype} column of
a \code{\link{gen_tibble}} object),
or a \code{\link{gen_tibble}}.}

\item{as_counts}{boolean defining whether the count of NAs (rather than the rate)
should be returned. It defaults to FALSE (i.e. rates are returned by default).}

\item{...}{currently unused.}
}
\value{
a vector of heterozygosities, one per individuals in the \code{\link{gen_tibble}}
}
\description{
Estimate missingness for each individual (i.e. the frequency of
missing genotypes in an individual).
}
