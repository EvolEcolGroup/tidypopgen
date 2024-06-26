% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop_fis.R
\name{pop_fis}
\alias{pop_fis}
\title{Compute population specific FIS}
\usage{
pop_fis(.x, include_global = FALSE, allele_sharing_mat = NULL)
}
\arguments{
\item{.x}{a grouped \code{\link{gen_tibble}} (as obtained by using \code{\link[dplyr:group_by]{dplyr::group_by()}})}

\item{include_global}{boolean determining whether, besides the population specific fis, a global
fis should be appended. Note that this will return a vector of n populations plus 1 (the global value)}

\item{allele_sharing_mat}{optional, the matrix of Allele Sharing returned by
\code{\link[=pairwise_allele_sharing]{pairwise_allele_sharing()}} with \code{as_matrix=TRUE}. As a number of statistics can be
derived from the Allele Sharing matrix,
it it sometimes more efficient to pre-compute this matrix.}
}
\value{
a vector of population specific fis (plus the global value if \code{include_global=TRUE})
}
\description{
This function computes population specific FIS (as computed by \code{\link[hierfstat:fs.dosage]{hierfstat::fis.dosage()}}).
}
