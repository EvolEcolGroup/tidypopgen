% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop_fst.R
\name{pop_fst}
\alias{pop_fst}
\title{Compute population specific Fst}
\usage{
pop_fst(.x, include_global = FALSE, allele_sharing_mat = NULL)
}
\arguments{
\item{.x}{a grouped \code{\link{gen_tibble}} (as obtained by using \code{\link[dplyr:group_by]{dplyr::group_by()}})}

\item{include_global}{boolean determining whether, besides the population
specific Fst, a global Fst should be appended. Note that this will return a
vector of n populations plus 1 (the global value)}

\item{allele_sharing_mat}{optional, the matrix of Allele Sharing returned by
\code{\link[=pairwise_allele_sharing]{pairwise_allele_sharing()}} with \code{as_matrix=TRUE}. As a number of
statistics can be derived from the Allele Sharing matrix,}
}
\value{
a vector of population specific Fst (plus the global value if
\code{include_global=TRUE})
}
\description{
This function computes population specific Fst, using the approach in Weir
and Goudet 2017 (as computed by \code{\link[hierfstat:fs.dosage]{hierfstat::fst.dosage()}}).
}
\references{
Weir, BS and Goudet J (2017) A Unified Characterization of
Population Structure and Relatedness. Genetics (2017) 206:2085
}
