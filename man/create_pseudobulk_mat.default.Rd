% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{create_pseudobulk_mat.default}
\alias{create_pseudobulk_mat.default}
\title{Create Pseudo-Bulk Matrix from scEpigenomics clusters & replicates}
\usage{
\method{create_pseudobulk_mat}{default}(
  object,
  by = "IDcluster",
  biological_replicate_col = NULL
)
}
\arguments{
\item{object}{A SingleCellExperiment object containing scEpigenomics dataset
with 'IDcluster' column.}

\item{by}{A character specifying the name of the metadata column referencing
the clusters.}

\item{biological_replicate_col}{Optional. A column of the
SingleCellExperiment object indicating the replicates or batches of the
dataset in order to take in account biological/technical noise. If NULL,
will create random layers of fake replicates.}
}
\value{
A pseudo-bulk matrix of cluster spread by replicates / batches /
fake replicates.
}
\description{
Create Pseudo-Bulk Matrix from scEpigenomics clusters & replicates
}
\examples{
if(requireNamespace("Seurat", quietly=TRUE)){
data("Seu", package = "IDclust")
mat <- create_pseudobulk_mat.default(object, by = "seurat_clusters")
}
}
