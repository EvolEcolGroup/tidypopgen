% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gt_save.R
\name{gt_get_file_names}
\alias{gt_get_file_names}
\title{Get the names of files storing the genotypes of a \code{gen_tibble}}
\usage{
gt_get_file_names(x)
}
\arguments{
\item{x}{a \code{\link{gen_tibble}}}
}
\value{
a character vector with the names and paths of the two files
}
\description{
A function to return the names of the files used to store data in a
\code{gen_tibble}. Specifically, this returns the .rds file storing the big
}
\examples{
example_gt <- load_example_gt("gen_tbl")

# To retrieve the names of and paths to the .bk and .rds files use:
gt_get_file_names(example_gt)

}
