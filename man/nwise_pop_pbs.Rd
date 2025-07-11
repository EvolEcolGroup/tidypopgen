% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nwise_pop_pbs.R
\name{nwise_pop_pbs}
\alias{nwise_pop_pbs}
\title{Compute the Population Branch Statistics for each combination of populations}
\usage{
nwise_pop_pbs(
  .x,
  type = c("tidy", "matrix"),
  fst_method = c("Hudson", "Nei87", "WC84"),
  return_fst = FALSE
)
}
\arguments{
\item{.x}{A grouped \code{gen_tibble}}

\item{type}{type of object to return. One of "tidy" or "matrix".
Default is "tidy".}

\item{fst_method}{the method to use for calculating Fst, one of 'Hudson',
'Nei87', and 'WC84'. See \code{\link[=pairwise_pop_fst]{pairwise_pop_fst()}} for details.}

\item{return_fst}{A logical value indicating whether to return the Fst values
along with the PBS values. Default is \code{FALSE}.}
}
\value{
Either a matrix with locus ID as rownames and the following columns:
\itemize{
\item \code{pbs_a.b.c}: the PBS value for population a given b & c (there
will be multiple such columns covering all 3 way combinations of
populations in the grouped \code{gen_tibble} object)
\item \code{pbsn1_a.b.c}: the normalized PBS value for population a given b & c.
\item \code{fst_a.b}: the Fst value for population a and b, if \code{return_fst} is TRUE
or a tidy tibble with the following columns:
\item \code{loci}: the locus ID
\item \code{stat_name}: the name of populations used in the pbs calculation
(e.g. "pbs_pop1.pop2.pop3"). If return_fst is TRUE, stat_name will also
include "fst" calculations in the same column (e.g. "fst_pop1.pop2").
\item \code{value}: the pbs value for the populations
}
}
\description{
The function computes the population branch statistics (PBS) for each
combination of populations at each locus. The PBS is a measure of the genetic
differentiation between one focal population and two reference populations,
and is used to identify outlier loci that may be under selection.
}
\examples{
example_gt <- load_example_gt()

# We can compute the PBS for all populations using "Hudson" method
example_gt \%>\%
  group_by(population) \%>\%
  nwise_pop_pbs(fst_method = "Hudson")
}
\references{
Yi X, et al. (2010) Sequencing of 50 human exomes reveals
adaptation to high altitude. Science 329: 75-78.
}
