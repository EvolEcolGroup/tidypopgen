% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_to_fbm_cpp.R
\name{vcf_to_fbm_cpp}
\alias{vcf_to_fbm_cpp}
\title{Convert vcf to FBM.}
\usage{
vcf_to_fbm_cpp(vcf_path, chunk_size = NULL, backingfile = NULL, quiet = FALSE)
}
\arguments{
\item{vcf_path}{the path to the vcf}

\item{chunk_size}{the chunk size to use on the vcf when loading the file}

\item{backingfile}{the name of the file to use as the backing file}
}
\value{
path to the resulting rds file as class bigSNP.
}
\description{
Convert a vcf file to a Filebacked Big Matrix (FBM) object.
This should work even for large vcf files that would not fit in memory.
}
\keyword{internal}
