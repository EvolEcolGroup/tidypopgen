% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_group_info.R
\name{annotate_group_info}
\alias{annotate_group_info}
\title{Annotate an admixture barplot with information about the populations}
\usage{
annotate_group_info(q_tbl, plt)
}
\arguments{
\item{q_tbl}{a tidied \code{q_matrix}
be in a block, all adjacent to each other}
}
\value{
modifier for a ggplot, added with the usual '+'
}
\description{
Add vertical lines to separate individuals from different populations,
and add optional x labels with their names
}
\keyword{internal}
