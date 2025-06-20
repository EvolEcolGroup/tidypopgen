% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyverse.R
\name{filter.gen_tbl}
\alias{filter.gen_tbl}
\title{Tidyverse methods for gt objects}
\usage{
\method{filter}{gen_tbl}(..., deparse.level = 1)
}
\arguments{
\item{...}{a gen_tibble and a data.frame or tibble}

\item{deparse.level}{an integer controlling the construction of column names.}
}
\value{
a \code{gen_tibble}
}
\description{
A filter method for \code{gen_tibble} objects
}
\examples{
test_gt <- load_example_gt("gen_tbl")
test_gt \%>\% filter(id \%in\% c("a", "c"))
test_gt <- load_example_gt("grouped_gen_tbl_sf")
test_gt \%>\% filter(id \%in\% c("a", "c"))
}
