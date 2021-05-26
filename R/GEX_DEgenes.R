#' Extracts the differentially expressed genes between two groups of cells. These groups are defined as cells having either of two entries (group1, group2) in the grouping.column of the input Seurat object metadata This function uses the FindMarkers function from the Seurat package. 
#' @param automate.GEX Output Seurat object from automate_GEX function that contained at least two distinct biological samples. The differential biological samples correspond to integer values in the order of the working directories initially supplied to the automate_GEX function.
#' @param min.pct The minimum percentage of cells expressing a gene in either of the two groups to be compared.
#' @param group1 either character or integer specifying the first group of cells that should be compared. (e.g. "s1" if sample_id is used as grouping.column)
#' @param group2 either character or integer specifying the first group of cells that should be compared. (e.g. "s2" if sample_id is used as grouping.column)
#' @param grouping.column Character. A column name of GEX@meta.data. In this column, group1 and group2 should be found. Defaults to "sample_id". Could also be set to "seurat_clusters" to generate DEGs between cells of 2 chosen clusters.
#' @param filter Character vector of initials of the genes to be filtered. Default is c("MT-", "RPL", "RPS"), which filters mitochondrial and ribosomal genes.
#' @param return.plot Logical specifying if a heatmap of the DEX genes is to be returned. If TRUE then @return is a list where the first element is a dataframe and the second a heatmap (see @return) 
#' @param logFC Logical specifying whether the genes will be displayed based on logFC (TRUE) or pvalue (FALSE).
#' @param up.genes Integer specifying the number of upregulated genes to be shown.
#' @param down.genes Integer specifying the number of downregulated genes to be shown. 
#' @param base The base with respect to which logarithms are computed. Default: 2
#' @param platypus.version. Function works with V2 and V3, no need to set this parameter
#' @return Returns a dataframe containing the output from the FindMarkers function, which contains information regarding the genes that are differentially regulated, statistics (p value and log fold change), and the percent of cells expressing the particular gene for both groups.
#' @export
#' @examples
#' \dontrun{
#' 
#' #Basic run between two samples
#' check_de_genes <- GEX_DEgenes(GEX= vgm[[2]],min.pct = .25,group1 = "s1",group2 = "s2")
#' 
#' #Getting DEGs between two seurat clusters
#' check_de_genes <- GEX_DEgenes(GEX= vgm[[2]],min.pct = .25, grouping_column = "seurat_clusters",group1 = "1",group2 = "4")
#' 
#' #Getting DEGs between two custom groups, possibly cellypes
#' check_de_genes <- GEX_DEgenes(GEX= vgm[[2]],min.pct = .25, grouping_column = "Column with cell type information",group1 = "T memory cells",group2 = "T effector cells")
#'}
GEX_DEgenes <- function(GEX, min.pct, group1, group2, grouping.column, filter, return.plot, logFC, up.genes, down.genes, base, platypus.version){
  
  if(missing(return.plot)) return.plot <- FALSE  
  if(missing(logFC)) logFC <- TRUE 
  if(missing(up.genes)) up.genes <- 15 
  if(missing(down.genes)) down.genes <- 15 
  if(missing(base)){base <- 2}
  platypus.version <- "does not matter"
  
  if (missing(filter)) filter <- c("MT-", "RPL", "RPS")
  
  if(missing(grouping.column)) grouping.column <- "sample_id"
  if(grouping.column %in% names(GEX@meta.data)){
    Seurat::Idents(GEX) <- GEX@meta.data[,c(grouping.column)]
  } else {
    stop("Please provide a valid GEX@meta.data column name as grouping.column")
  }

  cluster_markers <- Seurat::FindMarkers(GEX, min.pct = min.pct, ident.1 = as.character(group1), ident.2 = as.character(group2), base=base)
  colnames(cluster_markers)[2] <- "avg_logFC"
  cluster_markers$SYMBOL <- rownames(cluster_markers)
  
  exclude <- c()
  for (j in filter) {
    exclude <- c(exclude, stringr::str_which(rownames(cluster_markers), j))
  }
  cluster_markers <- cluster_markers[-exclude,]
  
  
  if (return.plot==TRUE) {
    if (logFC==TRUE) {
      ranks <- order(-cluster_markers$avg_logFC)
      cluster_markers <- cluster_markers[ranks,]
      heatmap_genes <- c(cluster_markers$SYMBOL[1:up.genes], cluster_markers$SYMBOL[(nrow(cluster_markers)-down.genes):nrow(cluster_markers)])
    }
    if (logFC==FALSE) {
      ranks <- order(cluster_markers$p_val_adj)
      cluster_markers <- cluster_markers[ranks,]
      heatmap_genes <- c(cluster_markers[which(cluster_markers$avg_logFC > 0),"SYMBOL"][1:up.genes], cluster_markers[which(cluster_markers$avg_logFC < 0),"SYMBOL"][1:down.genes])
    }
    cluster_markers_heatmap <- Seurat::DoHeatmap(GEX, features = heatmap_genes)
  }
  if (return.plot==FALSE) cluster_markers_heatmap <- NULL
  
  return(list(cluster_markers, cluster_markers_heatmap))
}
