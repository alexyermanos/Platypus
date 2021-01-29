
#' Extracts the differentially expressed genes between two samples. This function uses the FindMarkers function from the Seurat package. Further parameter control can be accomplished by calling the function directly on the output of automate_GEX and further extracting sample information from the "sample_id" component of the Seurat object.
#' @param automate.GEX Output Seurat object from automate_GEX function that contained at least two distinct biological samples. The differential biological samples correspond to integer values in the order of the working directories initially supplied to the automate_GEX function.
#' @param min.pct The minimum percentage of cells expressing a gene in either of the two groups to be compared.
#' @param sample1 either character or integer specifying the first sample that should be compared.
#' @param sample2 either character or integer specifying the first sample that should be compared.
#' @param by.group Logical specifying if groups should be used instead of samples. If TRUE, then the argument in sample1 and sample2 will correspond to cells found in the groups from sample1 or sample2.
#' @param filter Character vector of initials of the genes to be filtered. Default is c("MT-", "RPL", "RPS"), which filters mitochondrial and ribosomal genes.
#' @return Returns a dataframe containing the output from the FindMarkers function, which contains information regarding the genes that are differentially regulated, statistics (p value and log fold change), and the percent of cells expressing the particular gene for both groups.
#' @export
#' @examples
#' \dontrun{
#' check_de_genes2 <- GEX_DEgenes_persample(automate.GEX=automate.GEX.output[[i]],min.pct = .25,sample1 = "2",sample2 = "3")
#'}
GEX_DEgenes_persample <- function(automate.GEX, min.pct, sample1, sample2, by.group, filter){

    if(missing(by.group)) by.group <- FALSE
    if (missing(filter)) filter <- c("MT-", "RPL", "RPS")

    number_of_clusters <- length(unique(automate.GEX$seurat_clusters))
    if(by.group==TRUE) Seurat::Idents(automate.GEX) <- automate.GEX$group_id
    if(by.group==FALSE) Seurat::Idents(automate.GEX) <- automate.GEX$sample_id

    cluster_markers <- Seurat::FindMarkers(automate.GEX, min.pct = min.pct, ident.1 = as.character(sample1), ident.2 = as.character(sample2))
    cluster_markers$SYMBOL <- rownames(cluster_markers)
    exclude <- c()
    for (j in filter) {
        exclude <- c(exclude, !rownames(cluster_markers) %in% j)
    }
    cluster_markers <- cluster_markers[-exclude,]

    return(cluster_markers)
}
