#' Produces a heatmap displaying the expression of the top genes that define each cluster in the Seurat object. The output heatmap is derived from DoHeatmap from Seurat and thereby can be edited using typical ggplot interactions. The number of genes per cluster and the nunber of cells to display can be specified by the user. Either the log fold change or the p value can be used to select the top n genes.
#' @param automate_GEX.output Output Seurat object containing gene expression data from automate_GEX function that contained at least two distinct biological samples. The different biological samples correspond to integer values in the order of the working directories initially supplied to the automate_GEX function.
#' @param GEX_cluster_genes.output The output from the GEX_cluster_genes function - this should be a list with each list element corresponding to the genes, p values, logFC, pct expression for the genes differentially regulated for each cluster.
#' @param n.genes.per.cluster An integer value determining how many genes per cluster to display in the output heatmap. This number should be adjusted based on the number of clusters. Too many genes per cluster and clusters may cause a problem with the heatmap function in Seurat.
#' @param metric The metric that dictates which are the top n genes returned. Possible options are "p.value" (default) and "avg_logFC".
#' @param max.cell The max number of cells to display in the heatmap for each cluster, which corresponds to the number of columns. Default is set to 100 cells per cluster.
#' @return Returns a heatmap from the function DoHeatmap from the package Seurat, which is a ggplot object that can be modified or plotted. The number of genes is determined by the n.genes parameter and the number of cells per cluster is determined by the max.cell argument. This function gives a visual description of the top genes differentially expressed in each cluster.
#' @export
#' @examples
#' \dontrun{
#' cluster_defining_gene_heatmap <- GEX_cluster_genes_heatmap(automate_GEX.output=automate_GEX_output[[i]],GEX_cluster_genes.output=GEX_cluster_genes_output,n.genes.per.cluster=5,metric="p.value",max.cell=5)
#'}
GEX_cluster_genes_heatmap <- function(automate_GEX.output,GEX_cluster_genes.output,n.genes.per.cluster,metric,max.cell){
  if(missing(max.cell)) max.cell <- 100
  if(missing(n.genes.per.cluster)) n.genes.per.cluster <- 5
  if(missing(metric)) metric <- "p.value"

  holding_genes <- list()
  for(i in 1:length(GEX_cluster_genes.output)){
    if(metric=="p.value") holding_genes[[i]] <- rownames(GEX_cluster_genes.output[[i]][order(GEX_cluster_genes.output[[i]]$p_val_adj, decreasing = FALSE),])[1:n.genes.per.cluster]
    else if(metric=="avg_logFC") holding_genes[[i]] <- rownames(GEX_cluster_genes.output[[i]][order(abs(GEX_cluster_genes.output[[i]]$avg_logFC), decreasing = TRUE),])[1:n.genes.per.cluster]
  }
  ## Sample cells if too many
  sample_cells <- list()
  unique_clusters <- sort(unique(automate_GEX.output$seurat_clusters),decreasing = F)
  for(i in 1:length(unique_clusters)){
    if(length(which(automate_GEX.output$seurat_clusters==unique_clusters[i]))>max.cell){
      sample_cells[[i]] <- sample(which(automate_GEX.output$seurat_clusters==unique_clusters[i]),size = max.cell,replace = F)
    }
    else{
      sample_cells[[i]] <- which(automate_GEX.output$seurat_clusters==unique_clusters[i])
    }
  }
  output_heatmap <- Seurat::DoHeatmap(automate_GEX.output,features = unlist(unique(holding_genes)),cells = unlist(sample_cells))
  return(output_heatmap)
}
