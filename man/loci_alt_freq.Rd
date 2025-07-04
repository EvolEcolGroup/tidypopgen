% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loci_alt_freq.R
\name{loci_alt_freq}
\alias{loci_alt_freq}
\alias{loci_alt_freq.tbl_df}
\alias{loci_alt_freq.vctrs_bigSNP}
\alias{loci_alt_freq.grouped_df}
\alias{loci_maf}
\alias{loci_maf.tbl_df}
\alias{loci_maf.vctrs_bigSNP}
\alias{loci_maf.grouped_df}
\title{Estimate allele frequencies at each locus}
\usage{
loci_alt_freq(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores,
  block_size,
  type,
  ...
)

\method{loci_alt_freq}{tbl_df}(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)

\method{loci_alt_freq}{vctrs_bigSNP}(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
)

\method{loci_alt_freq}{grouped_df}(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  type = c("tidy", "list", "matrix"),
  ...
)

loci_maf(.x, .col = "genotypes", n_cores, block_size, type, ...)

\method{loci_maf}{tbl_df}(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)

\method{loci_maf}{vctrs_bigSNP}(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
)

\method{loci_maf}{grouped_df}(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  type = c("tidy", "list", "matrix"),
  ...
)
}
\arguments{
\item{.x}{a vector of class \code{vctrs_bigSNP} (usually the \code{genotypes} column of
a \code{\link{gen_tibble}} object), or a \code{\link{gen_tibble}}.}

\item{.col}{the column to be used when a tibble (or grouped tibble is passed
directly to the function). This defaults to "genotypes" and can only take
that value. There is no need for the user to set it, but it is included to
resolve certain tidyselect operations.}

\item{as_counts}{boolean defining whether the count of alternate and valid
(i.e. total number) alleles (rather than the frequencies) should be
returned. It defaults to FALSE (i.e. frequencies are returned by default).}

\item{n_cores}{number of cores to be used, it defaults to
\code{\link[bigstatsr:reexports]{bigstatsr::nb_cores()}}}

\item{block_size}{maximum number of loci read at once.}

\item{type}{type of object to return, if using grouped method. One of "tidy",
"list", or "matrix". Default is "tidy".}

\item{...}{other arguments passed to specific methods, currently unused.}
}
\value{
a vector of frequencies, one per locus, if \code{as_counts = FALSE};
else a matrix of two columns, the count of alternate alleles and the count
valid alleles (i.e. the sum of alternate and reference)
}
\description{
Allele frequencies can be estimates as minimum allele frequencies (MAF) with
\code{loci_maf()} or the frequency of the alternate allele (with
\code{loci_alt_freq()}). The latter are in line with the genotypes matrix (e.g. as
extracted by \code{\link[=show_loci]{show_loci()}}). Most users will be in interested in the MAF,
but the raw frequencies might be useful when computing aggregated statistics.
Both \code{loci_maf()} and \code{loci_alt_freq()} have efficient methods to support
grouped \code{gen_tibble} objects. These can return a tidied tibble, a list, or a
matrix.
}
\examples{
\dontshow{
data.table::setDTthreads(2)
RhpcBLASctl::blas_set_num_threads(2)
RhpcBLASctl::omp_set_num_threads(2)
}
example_gt <- load_example_gt("gen_tbl")

# For alternate allele frequency
example_gt \%>\% loci_alt_freq()

# For alternate allele frequency per locus per population
example_gt \%>\%
  group_by(population) \%>\%
  loci_alt_freq()
# alternatively, return a list of populations with their frequencies
example_gt \%>\%
  group_by(population) \%>\%
  loci_alt_freq(type = "list")
# or a matrix with populations in columns and loci in rows
example_gt \%>\%
  group_by(population) \%>\%
  loci_alt_freq(type = "matrix")
# or within reframe (not recommended, as it much less efficient
# than using it directly as shown above)
library(dplyr)
example_gt \%>\%
  group_by(population) \%>\%
  reframe(alt_freq = loci_alt_freq(genotypes))
# For MAF
example_gt \%>\% loci_maf()

# For minor allele frequency per locus per population
example_gt \%>\%
  group_by(population) \%>\%
  loci_maf()
# alternatively, return a list of populations with their frequencies
example_gt \%>\%
  group_by(population) \%>\%
  loci_maf(type = "list")
# or a matrix with populations in columns and loci in rows
example_gt \%>\%
  group_by(population) \%>\%
  loci_maf(type = "matrix")

}
