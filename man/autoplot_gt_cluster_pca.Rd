% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gt_cluster_pca.R
\name{autoplot.gt_cluster_pca}
\alias{autoplot.gt_cluster_pca}
\title{Autoplots for \code{gt_cluster_pca} objects}
\usage{
\method{autoplot}{gt_cluster_pca}(object, metric = c("BIC", "AIC", "WSS"), ...)
}
\arguments{
\item{object}{an object of class \code{gt_dapc}}

\item{metric}{the metric to plot on the y axis, one of 'BIC', 'AIC', or
'WSS' (with sum of squares)}

\item{...}{not currently used.}
}
\value{
a \code{ggplot2} object
}
\description{
For \code{gt_cluster_pca}, \code{autoplot} produces a plot of a metric of choice
('BIC', 'AIC' or 'WSS') against the number of clusters (\emph{k}). This plot is
can be used to infer the best value of \emph{k}, which corresponds to the smallest
value of the metric (the minimum in an 'elbow' shaped curve). In some cases,
there is not 'elbow' and the metric keeps decreasing with increasing \emph{k}; in
such cases, it is customary to choose the value of \emph{k} at which the decrease
in the metric reaches as plateau. For a programmatic way of choosing
\emph{k}, use \code{\link[=gt_cluster_pca_best_k]{gt_cluster_pca_best_k()}}.
}
\details{
\code{autoplot} produces simple plots to quickly inspect an object. They are not
customisable; we recommend that you use \code{ggplot2} to produce publication
ready plots.
}
\examples{
# Create a gen_tibble of lobster genotypes
bed_file <-
  system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
lobsters <- gen_tibble(bed_file,
  backingfile = tempfile("lobsters"),
  quiet = TRUE
)

# Remove monomorphic loci and impute
lobsters <- lobsters \%>\% select_loci_if(loci_maf(genotypes) > 0)
lobsters <- gt_impute_simple(lobsters, method = "mode")

# Create PCA object
pca <- gt_pca_partialSVD(lobsters)

# Run clustering on the first 10 PCs
cluster_pca <- gt_cluster_pca(
  x = pca,
  n_pca = 10,
  k_clusters = c(1, 5),
  method = "kmeans",
  n_iter = 1e5,
  n_start = 10,
  quiet = FALSE
)

# Autoplot BIC
autoplot(cluster_pca, metric = "BIC")

# # Autoplot AIC
autoplot(cluster_pca, metric = "AIC")

# # Autoplot WSS
autoplot(cluster_pca, metric = "WSS")
}
