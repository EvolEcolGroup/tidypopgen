% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_loci.R
\name{show_loci}
\alias{show_loci}
\alias{show_loci.tbl_df}
\alias{show_loci.vctrs_bigSNP}
\alias{show_loci<-}
\alias{show_loci<-.tbl_df}
\alias{show_loci<-.vctrs_bigSNP}
\title{Show the loci information of a \code{gen_tibble}}
\usage{
show_loci(.x, ...)

\method{show_loci}{tbl_df}(.x, ...)

\method{show_loci}{vctrs_bigSNP}(.x, ...)

show_loci(.x) <- value

\method{show_loci}{tbl_df}(.x) <- value

\method{show_loci}{vctrs_bigSNP}(.x) <- value
}
\arguments{
\item{.x}{a vector of class \code{vctrs_bigSNP} (usually the \code{genotype} column of
a \code{\link{gen_tibble}} object),
or a \code{\link{gen_tibble}}.}

\item{...}{currently unused.}

\item{value}{a data.frame or tibble of loci information to replace the
current one.}
}
\value{
a \code{\link[tibble:tibble]{tibble::tibble}} of information (see \code{\link{gen_tibble}} for details
on compulsory columns that will always be present)
}
\description{
Extract and set the information on loci from a  \code{gen_tibble}.
}
\examples{
example_gt <- load_example_gt("gen_tbl")

example_gt \%>\% show_loci()
}
