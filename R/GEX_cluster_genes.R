#' Extracts the differentially expressed genes between two samples. This function uses the FindMarkers function from the Seurat package. Further parameter control can be accomplished by calling the function directly on the output of automate_GEX or VDJ_GEX_matrix
#' @param GEX Output Seurat object of either automate_GEX for platypus.version v2 or of VDJ_GEX_matrix for platypus.version v3 (usually VDJ_GEX_matrix.output[[2]])
#' @param min.pct The minimum percentage of cells expressing a gene in either of the two groups to be compared. Default is 0.25
#' @param filter Character vector of initials of the genes to be filtered. Default is c("MT-", "RPL", "RPS"), which filters mitochondrial and ribosomal genes.
#' @param base The base with respect to which logarithms are computed. Default: 2
#' @param platypus.version is set automatically
#' @return Returns a dataframe containing the output from the FindMarkers function, which contains information regarding the genes that are differentially regulated, statistics (p value and log fold change), and the percent of cells expressing the particular gene. Ech element in the list corresponds to the clusters in numerical order. For example, the first element in the list output[[1]] corresponds to the genes deferentially expressed in cluster 0 in GEX
#' @export
#' @examples
#' \dontrun{
#' #Platypus version v2
#' GEX_cluster_genes(GEX =automate_GEX_output[[i]], min.pct = .25
#' , filter = c("MT-", "RPL", "RPS"))
#'
#' #Platypus version v3
#' GEX_cluster_genes(GEX = VDJ.GEX.matrix.output[[2]], min.pct = .25
#' , filter = c("MT-", "RPL", "RPS"))
#'}
GEX_cluster_genes <- function(GEX,
                              min.pct,
                              filter,
                              base,
                              platypus.version){

  platypus.version <- "does not matter"
  automate_GEX.output <- GEX
  GEX <- NULL

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

