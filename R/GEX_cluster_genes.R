#' Extracts the differentially expressed genes between two samples. This function uses the FindMarkers function from the Seurat package. Further parameter control can be accomplished by calling the function directly on the output of automate_GEX.
#' @param automate_GEX.output Output Seurat object containing gene expression data from automate_GEX function that contained at least two distinct biological samples. The different biological samples correspond to integer values in the order of the working directories initially supplied to the automate_GEX function.
#' @param min.pct The minimum percentage of cells expressing a gene in either of the two groups to be compared. Default is 0.25
#' @param filter Character vector of initials of the genes to be filtered. Default is c("MT-", "RPL", "RPS"), which filters mitochondrial and ribosomal genes.
#' @param base The base with respect to which logarithms are computed. Default: 2
#' @return Returns a dataframe containing the output from the FindMarkers function, which contains information regarding the genes that are differentially regulated, statistics (p value and log fold change), and the percent of cells expressing the particular gene. Ech element in the list corresponds to the clusters in numerical order. For example, the first element in the list output[[1]] corresponds to the genes differentially expressed in cluster 0 in the automate_GEX.output object.
#' @export
#' @examples
#' \dontrun{
#' genes_per_cluster <- GEX_cluster_genes(automate_GEX.output=automate_GEX_output[[i]], min.pct = .25, filter = c("MT-", "RPL", "RPS"))
#'}
GEX_cluster_genes <- function(automate_GEX.output, min.pct, filter, base){
  require(stringr)
  if(missing(min.pct)) min.pct <- 0.25
  if (missing(filter)) {filter <- c("MT-", "RPL", "RPS")}
  if(missing(base)){base <- 2}
  Seurat::Idents(automate_GEX.output) <- automate_GEX.output$seurat_clusters
  number_of_clusters <- length(unique(automate_GEX.output$seurat_clusters))
  cluster_markers <- list()
  for(i in 1:number_of_clusters){
    cluster_markers[[i]] <- Seurat::FindMarkers(automate_GEX.output, ident.1 = i-1, min.pct = min.pct, base=base)
    colnames(cluster_markers[[i]])[2] <- "avg_logFC"
    cluster_markers[[i]]$SYMBOL <- rownames(cluster_markers[[i]])
    cluster_markers[[i]]$cluster <- rep((i-1), nrow(cluster_markers[[i]]))
    exclude <- c()
    for (j in filter) {
      exclude <- c(exclude, stringr::str_which(rownames(cluster_markers[[i]]), j))
    }
    cluster_markers[[i]] <- cluster_markers[[i]][-exclude,]
  }
  return(cluster_markers)
}

