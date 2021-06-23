#' Organizes the top N genes that define each Seurat cluster and converts them into a single dataframe. This can be useful for obtaining insight into cluster-specific phenotypes.
#' @param GEX_cluster_genes.output The output from the GEX_cluster_genes function - this should be a list with each list element corresponding to the genes, p values, logFC, pct expression for the genes differentially regulated for each cluster.
#' @param n.genes The number of genes to be selected from each cluster. If n.genes is higher than the number of cells in a cluster then it is silently adjusted to be
#' @param by_FC Logical indicating if the top n genes are selected based on the logFC value instead of p value. Default is FALSE.
#' @param filter Character vector of initials of the genes to be filtered. Default is c("MT-", "RPL", "RPS"), which filters mitochondrial and ribosomal genes.
#' @return Returns a dataframe in which the top N genes defining each cluster based on differential expression are selected.
#' @importFrom dplyr %>%
#' @export
#' @examples
#' \dontrun{
#' GEX_topDE_genes_per_cluster(GEX_cluster_genes.output=list_of_genes_per_cluster
#' , n.genes=20, by_FC=FALSE, filter=c("MT-", "RPS", "RPL"))
#'}
GEX_topN_DE_genes_per_cluster <- function(GEX_cluster_genes.output, n.genes, by_FC, filter){

  abs_value <- NULL
  p_val <- NULL

  require(dplyr)
  require(stringr)
  if (missing(filter)) {filter <- c("MT-", "RPL", "RPS")}
  if (missing(by_FC)) {by_FC <- TRUE}

  output_list <- list()
  for(i in 1:length(GEX_cluster_genes.output)) {

    if("avg_log2FC" %in% names(GEX_cluster_genes.output[[i]])){
      names(GEX_cluster_genes.output[[i]])[which(names(GEX_cluster_genes.output[[i]]) == "avg_log2FC")] <- "avg_logFC"
    }

    temp.n.genes <- n.genes
    GEX_cluster_genes.output[[i]]$SYMBOL <- rownames(GEX_cluster_genes.output[[i]])
    GEX_cluster_genes.output[[i]]$cluster <- rep((i - 1),nrow(GEX_cluster_genes.output[[i]]))
    exclude <- c()
    for (j in filter) {
      exclude <- c(exclude, stringr::str_which(rownames(GEX_cluster_genes.output[[i]]),j))
    }
    if (length(exclude) > 0) {
      topN_filtered <- GEX_cluster_genes.output[[i]][-exclude,]
    }
    else {
      topN_filtered <- GEX_cluster_genes.output[[i]]
    }
    if (nrow(topN_filtered) < n.genes) {
      temp.n.genes <- nrow(topN_filtered)
    }
    if (by_FC) {
      topN_filtered <- topN_filtered %>% dplyr::mutate(abs_value = abs(topN_filtered$avg_logFC))
      output_list[[i]] <- dplyr::select(dplyr::slice_max(topN_filtered,n = temp.n.genes, abs_value),!abs_value)
    }
    else {
      output_list[[i]] <- dplyr::slice_min(topN_filtered,
                                           n = temp.n.genes, p_val)
    }
  }
  output_unlist <- do.call("rbind", output_list)
  return(output_unlist)
}
