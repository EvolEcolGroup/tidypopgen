% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gt_extract_f2.R
\name{gt_extract_f2}
\alias{gt_extract_f2}
\title{Compute and store blocked f2 statistics for ADMIXTOOLS 2}
\usage{
gt_extract_f2(
  .x,
  outdir = NULL,
  blgsize = 0.05,
  maxmem = 8000,
  maxmiss = 0,
  minmaf = 0,
  maxmaf = 0.5,
  minac2 = FALSE,
  outpop = NULL,
  outpop_scale = TRUE,
  transitions = TRUE,
  transversions = TRUE,
  overwrite = FALSE,
  adjust_pseudohaploid = NULL,
  fst = TRUE,
  afprod = TRUE,
  poly_only = c("f2"),
  apply_corr = TRUE,
  n_cores = 1,
  quiet = FALSE
)
}
\arguments{
\item{.x}{a \code{\link{gen_tibble}}}

\item{outdir}{Directory where data will be stored.}

\item{blgsize}{SNP block size in Morgan. Default is 0.05 (5 cM). If \code{blgsize}
is 100 or greater, if will be interpreted as base pair distance rather than
centimorgan distance.}

\item{maxmem}{Maximum amount of memory to be used. If the required amount of
memory exceeds \code{maxmem}, allele frequency data will be split into blocks,
and the computation will be performed separately on each block pair. This
doesn't put a precise cap on the amount of memory used (it used to at some
point). Set this parameter to lower values if you run out of memory while
running this function. Set it to higher values if this function is too slow
and you have lots of memory.}

\item{maxmiss}{Discard SNPs which are missing in a fraction of populations
higher than \code{maxmiss}}

\item{minmaf}{Discard SNPs with minor allele frequency less than \code{minmaf}}

\item{maxmaf}{Discard SNPs with minor allele frequency greater than than
\code{maxmaf}}

\item{minac2}{Discard SNPs with allele count lower than 2 in any population
(default \code{FALSE}). This option should be set to \code{TRUE} when computing
f3-statistics where one population consists mostly of pseudohaploid
samples. Otherwise heterozygosity estimates and thus f3-estimates can be
biased. \code{minac2 == 2} will discard SNPs with allele count lower than 2 in
any non-singleton population (this option is experimental and is based on
the hypothesis that using SNPs with allele count lower than 2 only leads to
biases in non-singleton populations). Note that, While the \code{minac2} option
discards SNPs with allele count lower than 2 in any population, the
\code{qp3pop} function will only discard SNPs with allele count lower than
2 in the first (target) population (when the first argument is the prefix
of a genotype file; i.e. it is applied directly to a genotype file, not via
precomputing f2 from a \code{\link{gen_tibble}}).}

\item{outpop}{Keep only SNPs which are heterozygous in this population}

\item{outpop_scale}{Scale f2-statistics by the inverse \code{outpop}
heterozygosity (\code{1/(p*(1-p))}). Providing \code{outpop} and setting
\code{outpop_scale} to \code{TRUE} will give the same results as the original
\emph{qpGraph} when the \code{outpop} parameter has been set, but it has the
disadvantage of treating one population different from the others. This may
limit the use of these f2-statistics for other models.}

\item{transitions}{Set this to \code{FALSE} to exclude transition SNPs}

\item{transversions}{Set this to \code{FALSE} to exclude transversion SNPs}

\item{overwrite}{Overwrite existing files in \code{outdir}}

\item{adjust_pseudohaploid}{Genotypes of pseudohaploid samples are usually
coded as \code{0} or \code{2}, even though only one allele is observed.
\code{adjust_pseudohaploid} ensures that the observed allele count increases
only by \code{1} for each pseudohaploid sample. If \code{TRUE} (default), samples
that don't have any genotypes coded as \code{1} among the first 1000 SNPs are
automatically identified as pseudohaploid. This leads to slightly more
accurate estimates of f-statistics. Setting this parameter to \code{FALSE}
treats all samples as diploid and is equivalent to the \emph{ADMIXTOOLS} \code{ inbreed: NO} option. Setting \code{adjust_pseudohaploid} to an integer \code{n} will
check the first \code{n} SNPs instead of the first 1000 SNPs. NOW DEPRECATED,
set the ploidy of the \code{gen_tibble} with \code{\link[=gt_pseudohaploid]{gt_pseudohaploid()}}.}

\item{fst}{Write files with pairwise FST for every population pair. Setting
this to FALSE can make \code{extract_f2} faster and will require less memory.}

\item{afprod}{Write files with allele frequency products for every
population pair. Setting this to FALSE can make \code{extract_f2} faster and
will require less memory.}

\item{poly_only}{Specify whether SNPs with identical allele frequencies in
every population should be discarded (\code{poly_only = TRUE}), or whether they
should be used (\code{poly_only = FALSE}). By default (\code{poly_only = c("f2")}),
these SNPs will be used to compute FST and allele frequency products, but
not to compute f2 (this is the default option in the original ADMIXTOOLS).}

\item{apply_corr}{Apply small-sample-size correction when computing
f2-statistics (default \code{TRUE})}

\item{n_cores}{Parallelize computation across \code{n_cores} cores.}

\item{quiet}{Suppress printing of progress updates}
}
\value{
SNP metadata (invisibly)
}
\description{
This function prepares data for various \emph{ADMIXTOOLS 2} functions from the
package \emph{ADMIXTOOLS 2}. It takes a \code{\link{gen_tibble}}, computes allele
frequencies and blocked f2-statistics, and writes the results to \code{outdir}. It
is equivalent to \code{admixtools::extract_f2()}.
}
\examples{
# run the example only if we have the package installed
\donttest{
if (requireNamespace("admixtools", quietly = TRUE)) {
  bed_file <-
    system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
  lobsters <- gen_tibble(bed_file,
    backingfile = tempfile("lobsters"),
    quiet = TRUE
  )
  lobsters <- lobsters \%>\% group_by(population)
  f2_path <- tempfile()
  gt_extract_f2(lobsters, outdir = f2_path, quiet = TRUE)
  admixtools::f2_from_precomp(f2_path, verbose = FALSE)
}
}
}
