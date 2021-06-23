#' Extracts the differentially expressed genes between two samples. This function uses the FindMarkers function from the Seurat package. Further parameter control can be accomplished by calling the function directly on the output of automate_GEX and further extracting sample information from the "sample_id" component of the Seurat object.
#' @param automate.GEX Output Seurat object from automate_GEX function that contained at least two distinct biological samples. The differential biological samples correspond to integer values in the order of the working directories initially supplied to the automate_GEX function.
#' @param min.pct The minimum percentage of cells expressing a gene in either of the two groups to be compared.
#' @param sample1 either character or integer specifying the first sample that should be compared.
#' @param sample2 either character or integer specifying the first sample that should be compared.
#' @param by.group Logical specifying if groups should be used instead of samples. If TRUE, then the argument in sample1 and sample2 will correspond to cells found in the groups from sample1 or sample2.
#' @param filter Character vector of initials of the genes to be filtered. Default is c("MT-", "RPL", "RPS"), which filters mitochondrial and ribosomal genes.
#' @param return.plot Logical specifying if a heatmap of the DEX genes is to be returned. If TRUE then @return is a list where the first element is a dataframe and the second a heatmap (see @return)
#' @param logFC Logical specifying whether the genes will be displayed based on logFC (TRUE) or pvalue (FALSE).
#' @param up.genes Integer specifying the number of upregulated genes to be shown.
#' @param down.genes Integer specifying the number of downregulated genes to be shown.
#' @param base The base with respect to which logarithms are computed. Default: 2
#' @return Returns a dataframe containing the output from the FindMarkers function, which contains information regarding the genes that are differentially regulated, statistics (p value and log fold change), and the percent of cells expressing the particular gene for both groups.
#' @export
#' @examples
#' \dontrun{
#' check_de_genes2 <- GEX_DEgenes_persample(automate.GEX=automate.GEX.output[[i]]
#' ,min.pct = .25,sample1 = "1",sample2 = "2")
#'}
GEX_DEgenes_persample <- function(automate.GEX, min.pct, sample1, sample2, by.group, filter, return.plot, logFC, up.genes, down.genes, base){

    if(missing(return.plot)) return.plot <- FALSE
    if(missing(logFC)) logFC <- TRUE
    if(missing(up.genes)) up.genes <- 15
    if(missing(down.genes)) down.genes <- 15
    if(missing(base)){base <- 2}

    if(missing(by.group)) by.group <- FALSE
    if (missing(filter)) filter <- c("MT-", "RPL", "RPS")

    number_of_clusters <- length(unique(automate.GEX$seurat_clusters))
    if(by.group==TRUE) Seurat::Idents(automate.GEX) <- automate.GEX$group_id
    if(by.group==FALSE) Seurat::Idents(automate.GEX) <- automate.GEX$sample_id

    cluster_markers <- Seurat::FindMarkers(automate.GEX, min.pct = min.pct, ident.1 = as.character(sample1), ident.2 = as.character(sample2), base=base)
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
      cluster_markers_heatmap <- Seurat::DoHeatmap(automate.GEX, features = heatmap_genes)
    }
    if (return.plot==FALSE) cluster_markers_heatmap <- NULL

    return(list(cluster_markers, cluster_markers_heatmap))
}
