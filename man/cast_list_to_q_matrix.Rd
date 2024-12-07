% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/q_matrix.R
\name{cast_list_to_q_matrix}
\alias{cast_list_to_q_matrix}
\title{Creates a structured list of Q matrices and their K indices}
\usage{
cast_list_to_q_matrix(matrix_list)
}
\arguments{
\item{matrix_list}{a flat list of \code{q_matrix} objects}
}
\value{
a \code{q_matrix_list} object containing the Q matrices and K indices
}
\description{
This function takes two or more \code{q_matrix} objects in a flat list format and
returns them in a structured list \code{q_matrix_list} where the first element
k_indices contains the index number for each Q matrix grouped by K, and the
second element q_matrices contains the Q matrices themselves.
}
\keyword{internal}
