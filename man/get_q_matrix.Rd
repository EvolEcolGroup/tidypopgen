% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/q_matrix.R
\name{get_q_matrix}
\alias{get_q_matrix}
\title{Return a single Q matrix from a \code{gt_admix} object}
\usage{
get_q_matrix(x, ..., k, run)
}
\arguments{
\item{x}{A \code{gt_admix} object containing multiple Q matrices}

\item{...}{Not used}

\item{k}{The k value of the desired Q matrix}

\item{run}{The run number of the desired Q matrix}
}
\value{
A single Q matrix from the \code{gt_admix} object
}
\description{
This function retrieves a single Q matrix from a \code{gt_admix} object
based on the specified k value and run number.
}
