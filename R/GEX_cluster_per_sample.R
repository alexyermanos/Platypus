#' Plots the cluster membership for each of the distinct samples in the Seurat object from the automate_GEX function. The distinct samples are determined by "sample_id" field in the Seurat object.
#' @param automate_GEX.output Output Seurat object containing gene expression data from automate_GEX function that contained at least two distinct biological samples. The different biological samples correspond to integer values in the order of the working directories initially supplied to the automate_GEX function.
#' @return Returns a ggplot in which the values on the x axis correspond to each cluster found in the Seurat object. The y axis corresponds to the percentage of cells found in each cluster. The bar and color corresponds to the distinct sample_id.
#' @export
#' @examples
#' \dontrun{
#' cluster.distribution.per.sample <- GEX_cluster_per_sample(automate_GEX.output=automate_GEX_out[[i]])
#'}
GEX_cluster_per_sample <- function(automate_GEX.output){
  unique_samples <- sort(unique(automate_GEX.output$sample_id),decreasing = F)
  unique_clusters <- sort(unique(automate_GEX.output$seurat_clusters),decreasing = F)
  cells_per_cluster_per_sample <- list()
  for(i in 1:length(unique_samples)){
    cells_per_cluster_per_sample[[i]] <- list()
    for(j in 1:length(unique_clusters)){
      cells_per_cluster_per_sample[[i]][[j]] <- length(which(automate_GEX.output$sample_id==unique_samples[i] & automate_GEX.output$seurat_clusters==unique_clusters[j]))/length(which(automate_GEX.output$sample_id==unique_samples[i]))
    }
  }
  melting <- reshape2::melt(cells_per_cluster_per_sample)
  melting$L1 <- as.character(melting$L1)
  output.plot <- ggplot2::ggplot(melting, ggplot2::aes(fill = L1, y=value, x=L2,group=L1)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black",position = "dodge") + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Fraction of cells") + ggplot2::xlab("Cluster membership")
  return(output.plot)
}
