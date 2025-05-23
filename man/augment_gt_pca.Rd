% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gt_pca_tidiers.R
\name{augment_gt_pca}
\alias{augment_gt_pca}
\alias{augment.gt_pca}
\title{Augment data with information from a gt_pca object}
\usage{
\method{augment}{gt_pca}(x, data = NULL, k = NULL, ...)
}
\arguments{
\item{x}{A \code{gt_pca} object returned by one of the \verb{gt_pca_*} functions.}

\item{data}{the \code{gen_tibble} used to run the PCA.}

\item{k}{the number of components to add}

\item{...}{Not used. Needed to match generic signature only.}
}
\value{
A  \link{gen_tibble} containing the original data along with additional
columns containing each observation's projection into PCA space.
}
\description{
Augment for \code{gt_pca} accepts a model object and a dataset and adds scores to
each observation in the dataset. Scores for each component are stored in a
separate column, which is given name with the pattern ".fittedPC1",
".fittedPC2", etc. For consistency with \link[broom:augment.prcomp]{broom::augment.prcomp}, a column
".rownames" is also returned; it is a copy of 'id', but it ensures that any
scripts written for data augmented with \link[broom:augment.prcomp]{broom::augment.prcomp} will work out
of the box (this is especially helpful when adapting plotting scripts).
}
\seealso{
\code{\link[=gt_pca_autoSVD]{gt_pca_autoSVD()}} \link{gt_pca_tidiers}
}
