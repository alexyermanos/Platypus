#'Heatmap of cluster defining genes
#'
#'@description Produces a heatmap displaying the expression of the top genes that define each cluster in the Seurat object. The output heatmap is derived from DoHeatmap from Seurat and thereby can be edited using typical ggplot interactions. The number of genes per cluster and the nunber of cells to display can be specified by the user. Either the log fold change or the p value can be used to select the top n genes.
#' @param GEX Output Seurat object of either automate_GEX for platypus.version v2 or of VDJ_GEX_matrix for platypus.version v3 (usually VDJ_GEX_matrix.output[[2]])
#' @param GEX_cluster_genes.output The output from the GEX_cluster_genes function - this should be a list with each list element corresponding to the genes, p values, logFC, pct expression for the genes deferentially regulated for each cluster.
#' @param n.genes.per.cluster An integer value determining how many genes per cluster to display in the output heatmap. This number should be adjusted based on the number of clusters. Too many genes per cluster and clusters may cause a problem with the heatmap function in Seurat.
#' @param metric The metric that dictates which are the top n genes returned. Possible options are "p.value" (default), "avg_logFC", "top_logFC", "bottom_logFC". "top_logFC" returns the top expressed genes for each cluster, whereas "bottom_logFC" returns the least expressed genes per cluster-both by log fold change.
#' @param max.cell The max number of cells to display in the heatmap for each cluster, which corresponds to the number of columns. Default is set to 100 cells per cluster.
#' @param group.colors Optional character vector. Array of colors with the same length as GEX_cluster_genes.output to color bars above the heatmap. Defaults to rainbow palette
#' @param slot Seurat object slot from which to plot gene expression data.
#' @param platypus.version is set automatically
#' @return Returns a heatmap from the function DoHeatmap from the package Seurat, which is a ggplot object that can be modified or plotted. The number of genes is determined by the n.genes parameter and the number of cells per cluster is determined by the max.cell argument. This function gives a visual description of the top genes differentially expressed in each cluster.
#' @export
#' @examples
#'\donttest{
#'try({
#' GEX_cluster_genes_output <- GEX_cluster_genes(GEX =
#' subset(Platypus::small_vgm[[2]],
#' seurat_clusters %in% c(0,1)), min.pct = .25
#' , filter = c("MT-", "RPL", "RPS"))
#'
#' cluster_defining_gene_heatmap <- GEX_cluster_genes_heatmap(GEX =
#' Platypus::small_vgm[[2]]
#' ,GEX_cluster_genes.output=GEX_cluster_genes_output
#' ,n.genes.per.cluster=5,metric="p.value",max.cell=5)
#' })
#'}

GEX_cluster_genes_heatmap <- function(GEX,
                                      GEX_cluster_genes.output,
                                      n.genes.per.cluster,
                                      metric,
                                      max.cell,
                                      group.colors,
                                      slot,
                                      platypus.version){

  platypus.version <- "does not matter"
  if(missing(max.cell)) max.cell <- 100
  if(missing(n.genes.per.cluster)) n.genes.per.cluster <- 5
  if(missing(metric)) metric <- "p.value"
  if(missing(group.colors)) group.colors <- grDevices::rainbow(length(GEX_cluster_genes.output))
  if(missing(slot)) slot <- "counts"

  #rename in case of naming change
  if(any("avg_log2FC" %in% names(GEX_cluster_genes.output[[1]]))){
    fc_col <- which(names(GEX_cluster_genes.output[[1]]) == "avg_log2FC")
    for(i in 1:length(GEX_cluster_genes.output)){
      names(GEX_cluster_genes.output[[i]])[fc_col] <- "avg_logFC"
    }
  }

  holding_genes <- list()
  for(i in 1:length(GEX_cluster_genes.output)){
    if(metric=="p.value") holding_genes[[i]] <- rownames(GEX_cluster_genes.output[[i]][order(GEX_cluster_genes.output[[i]]$p_val_adj, decreasing = FALSE),])[1:n.genes.per.cluster]
    else if(metric=="avg_logFC") holding_genes[[i]] <- rownames(GEX_cluster_genes.output[[i]][order(abs(GEX_cluster_genes.output[[i]]$avg_logFC), decreasing = TRUE),])[1:n.genes.per.cluster]
    else if(metric=="top_logFC") holding_genes[[i]] <- rownames(GEX_cluster_genes.output[[i]][order((GEX_cluster_genes.output[[i]]$avg_logFC), decreasing = TRUE),])[1:n.genes.per.cluster]
    else if(metric=="bottom_logFC") holding_genes[[i]] <- rownames(GEX_cluster_genes.output[[i]][order((GEX_cluster_genes.output[[i]]$avg_logFC), decreasing = FALSE),])[1:n.genes.per.cluster]
  }
  ## Sample cells if too many
  sample_cells <- list()
  unique_clusters <- sort(unique(GEX$seurat_clusters),decreasing = FALSE)
  for(i in 1:length(unique_clusters)){
    if(length(which(GEX$seurat_clusters==unique_clusters[i]))>max.cell){
      sample_cells[[i]] <- sample(which(GEX$seurat_clusters==unique_clusters[i]),size = max.cell,replace = FALSE)
    }
    else{
      sample_cells[[i]] <- which(GEX$seurat_clusters==unique_clusters[i])
    }
  }
  output_heatmap <- Seurat::DoHeatmap(GEX,features = unlist(unique(holding_genes)),cells = unlist(sample_cells), group.colors = group.colors, slot = slot)
  return(output_heatmap)
}
