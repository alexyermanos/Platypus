% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GEX_pairwise_degs.R
\name{GEX_pairwise_DEGs}
\alias{GEX_pairwise_DEGs}
\title{Wrapper for calculating pairwise differentially expressed genes}
\usage{
GEX_pairwise_DEGs(
  GEX,
  group.by,
  min.pct,
  RP.MT.filter,
  label.n.top.genes,
  genes.to.label,
  save.plot,
  save.csv
)
}
\arguments{
\item{GEX}{Output Seurat object of the VDJ_GEX_matrix function (VDJ_GEX_matrix.output[[2]])}

\item{group.by}{Character. Defaults to "seurat_clusters" Column name of GEX@meta.data to use for pairwise comparisons. More than 20 groups are discuraged.}

\item{min.pct}{Numeric. Defaults to 0.25 passed to Seurat::FindMarkers}

\item{RP.MT.filter}{Boolean. Defaults to True. If True, mitochondrial and ribosomal genes are filtered out from the output of Seurat::FindMarkers}

\item{label.n.top.genes}{Integer. Defaults to 50. Defines how many genes are labelled via geom_text_repel. Genes are ordered by adjusted p value and the first label.n.genes are labelled}

\item{genes.to.label}{Character vector. Defaults to "none". Vector of gene names to plot independently of their p value. Can be used in combination with label.n.genes.}

\item{save.plot}{Boolean. Defaults to False. Whether to save plots as appropriately named .png files}

\item{save.csv}{Boolean. Defaults to False. Whether to save deg tables as appropriately named .csv files}
}
\value{
A nested list with out[[i]][[1]] being ggplot volcano plots and out[[i]][[2]] being source DEG dataframes.
}
\description{
Produces and saves a list of volcano plots with each showing differentially expressed genes between pairs groups. If e.g. seurat_clusters used as group.by, a plot will be generated for every pairwise comparison of clusters. For large numbers of this may take longer to run. Only available for platypus v3
}
\examples{
\dontrun{
GEX_pairwise_DEGs(GEX = Platypus::small_vgm[[2]],group.by = "sample_id"
,min.pct = 0.25,RP.MT.filter = TRUE,label.n.top.genes = 2,genes.to.label = c("CD24A")
,save.plot = FALSE, save.csv = FALSE)
}
}
