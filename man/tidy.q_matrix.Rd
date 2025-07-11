% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/q_matrix.R
\name{tidy.q_matrix}
\alias{tidy.q_matrix}
\title{Tidy a Q matrix}
\usage{
\method{tidy}{q_matrix}(x, data, ...)
}
\arguments{
\item{x}{A Q matrix object (as returned by \code{\link{q_matrix}}).}

\item{data}{An associated tibble (e.g. a \code{\link{gen_tibble}}), with the
individuals in the same order as the data used to generate the Q matrix}

\item{...}{not currently used}
}
\value{
A tidied tibble containing columns:
\item{\code{row}}{ID of the original observation (i.e. rowname from original
data).}
\item{\code{Q}}{Integer indicating a Q component.} \item{\code{value}}{The proportion
for that particular Q value.}
}
\description{
Takes a \code{q_matrix} object, which is a matrix, and returns a tidied tibble.
}
\examples{
# run the example only if we have the package installed
if (requireNamespace("LEA", quietly = TRUE)) {
  example_gt <- load_example_gt("gen_tbl")

  # Create a gt_admix object
  admix_obj <- example_gt \%>\% gt_snmf(k = 1:3, project = "force")

  # Extract a Q matrix
  q_mat_k3 <- get_q_matrix(admix_obj, k = 3, run = 1)

  tidy(q_mat_k3, data = example_gt)
}
}
