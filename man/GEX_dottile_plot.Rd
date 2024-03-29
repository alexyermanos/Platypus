% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GEX_dottile_plot.R
\name{GEX_dottile_plot}
\alias{GEX_dottile_plot}
\title{GEX Dottile plots}
\usage{
GEX_dottile_plot(GEX, genes, group.by, threshold.to.plot, platypus.version)
}
\arguments{
\item{GEX}{GEX seurat object generated with VDJ_GEX_matrix}

\item{genes}{Character vector. Genes of those in rownames(GEX) to plot. Can be any number, but more then 30 is discuraged because of cluttering}

\item{group.by}{Character. Name of a column in GEX@meta.data to split the plot by. If set to \"none\", a plot with a single column will be produced.}

\item{threshold.to.plot}{Integer 1-100. \% of cells which must be expressing the feature to plot a point. If below, the field will be left empty}

\item{platypus.version}{This is coded for \"v3\" only, but in practice any Seurat Object can be fed in}
}
\value{
Returns a ggplot object were the dot size indicates the percentage of expressing cells and the dot color indicates the expression level.
}
\description{
Outputs a dotplot for gene expression, where the color of each dot is scaled by the gene expression level and the size is scaled by the \% of cells positive for the gene
}
\examples{
#To return a plot detailing the expression of common genes by seurat cluster
GEX_dottile_plot(GEX = Platypus::small_vgm[[2]], genes = c("CD19","CD83"),
group.by = "seurat_clusters", threshold.to.plot = 5)

}
