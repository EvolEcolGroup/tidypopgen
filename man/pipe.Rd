% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\arguments{
\item{lhs}{A value or the magrittr placeholder.}

\item{rhs}{A function call using the magrittr semantics.}
}
\value{
The result of calling \code{rhs(lhs)}.
}
\description{
See `magrittr::pipe \\%>\\% for details.
}
\examples{
example_gt <- load_example_gt("gen_tbl")
example_gt \%>\% count_loci()
}
\keyword{internal}
