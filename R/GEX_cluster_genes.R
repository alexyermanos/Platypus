#' Extracts the differentially expressed genes for each cluster as determined by Seurat. This function uses the FindMarkers function from the Seurat package. Further parameter control can be accomplished by calling the function directly on the output of automate_GEX.
#' @param automate_GEX.output Output Seurat object containing gene expression data from automate_GEX function that contained at least two distinct biological samples. The different biological samples correspond to integer values in the order of the working directories initially supplied to the automate_GEX function.
#' @param min.pct The minimum percentage of cells expressing a gene in either of the two groups to be compared. Default is 0.25
#' @return Returns a dataframe containing the output from the FindMarkers function, which contains information regarding the genes that are differentially regulated, statistics (p value and log fold change), and the percent of cells expressing the particular gene. Ech element in the list corresponds to the clusters in numerical order. For example, the first element in the list output[[1]] corresponds to the genes differentially expressed in cluster 0 in the automate_GEX.output object.
#' @export
#' @examples
#' \dontrun{
#' genes_per_cluster <- GEX_cluster_genes(automate_GEX.output=automate_GEX_output[[i]],min.pct = .25)
#'}
GEX_cluster_genes <- function(automate_GEX.output,min.pct){
  if(missing(min.pct)) min.pct <- 0.25
  Seurat::Idents(automate_GEX.output) <- automate_GEX.output$seurat_clusters
  number_of_clusters <- length(unique(automate_GEX.output$seurat_clusters))
  cluster_markers <- list()
  for(i in 1:number_of_clusters){
    cluster_markers[[i]] <- Seurat::FindMarkers(automate_GEX.output, ident.1 = i-1, min.pct = min.pct)
  }
  return(cluster_markers)
}
