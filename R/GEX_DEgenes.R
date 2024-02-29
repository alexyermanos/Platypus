#'Wrapper for differential gene expression analysis and plotting
#'
#'@description Extracts the differentially expressed genes between two groups of cells. These groups are defined as cells having either of two entries (group1, group2) in the grouping.column of the input Seurat object metadata This function uses the FindMarkers function from the Seurat package.
#' @param GEX Output Seurat object from automate_GEX or VDJ_GEX_matrix_function (VDJ_GEX_matrix.output[[2]]) function that contained at least two distinct biological groups.
#' @param FindMarkers.out OPTIONAL: the output of the FindMarkers function. This skips the DEG calculation step and outputs desired plots. All plotting parameters function as normal. Grouping parameters and min.pct are ignored.
#' @param grouping.column Character. A column name of GEX@meta.data. In this column, group1 and group2 should be found. Defaults to "sample_id". Could also be set to "seurat_clusters" to generate DEGs between cells of 2 chosen clusters.
#' @param group1 either character or integer specifying the first group of cells that should be compared. (e.g. "s1" if sample_id is used as grouping.column)
#' @param group2 either character or integer specifying the first group of cells that should be compared. (e.g. "s2" if sample_id is used as grouping.column)
#' @param filter Character vector of initials of the genes to be filtered. Default is c("MT-", "RPL", "RPS"), which filters mitochondrial and ribosomal genes.
#' @param min.pct The minimum percentage of cells expressing a gene in either of the two groups to be compared.
#' @param base The base with respect to which logarithms are computed. Default: 2
#' @param logFC Logical specifying whether the genes will be displayed based on logFC (TRUE) or pvalue (FALSE).
#' @param color.p.threshold numeric specifying the adjusted p-value threshold for geom_points to be colored. Default is set to 0.01.
#' @param color.log.threshold numeric specifying the absolute logFC threshold for geom_points to be colored. Default is set to 0.25.
#' @param color.by.threshold Boolean. Set to TRUE to color by color.p.threshold and color.log.threshold. Set to FALSE for a continuous color scale by fold change.
#' @param return.plot Character specifying if a "heatmap", "heatmap" or a "volcano" or "none" is to be returned. If not "none" then @return is a list where the first element is a dataframe and the second a plot (see @return). Defaults to none
#' @param up.genes FOR HEATMAP Integer specifying the number of upregulated genes to be shown.
#' @param down.genes FOR HEATMAP Integer specifying the number of downregulated genes to be shown.
#' @param genes.to.label FOR VOLCANO Character vector of genes to label irregardless of their p value.
#' @param label.n.top.genes FOR VOLCANO Interger. How many top genes to label either by Fold change (if logFC ==TRUE) or by p.value (if logFC == FALSE). More than 50 are not recommended. Also works in conjunction with genes.to.label
#' @param platypus.version Function works with V2 and V3, no need to set this parameter
#' @param size.top.colorbar Integer. Size of the top colorbar for heatmap plot.
#' @return Returns a dataframe containing the output from the FindMarkers function, which contains information regarding the genes that are differentially regulated, statistics (p value and log fold change), and the percent of cells expressing the particular gene for both groups.
#' @export
#' @examples
#' DEGs <- GEX_DEgenes(GEX = Platypus::small_vgm[[2]],min.pct = .25,
#' group1 = "s1",group2 = "s2", return.plot = "volcano")


GEX_DEgenes <- function(GEX,
                           FindMarkers.out,
                           grouping.column,
                           group1,
                           group2,
                           min.pct,
                           filter,
                           return.plot,
                           logFC,
                           color.p.threshold,
                           color.log.threshold,
                           color.by.threshold,
                           up.genes,
                           down.genes,
                           base,
                           label.n.top.genes,
                           genes.to.label,
                           platypus.version,
                           size.top.colorbar){
  #Variation of the function so you can plot the title
  SYMBOL <- NULL
  avg_logFC <- NULL
  p_val_adj <- NULL
  color_threshold <- NULL

  if(missing(FindMarkers.out)) FindMarkers.out <- "none"
  if(missing(return.plot)) return.plot <- "none"
  if(return.plot == TRUE){
    message("Please set return.plot to either 'volcano', 'heatmap' or 'none'")
    message("Setting return.plot to 'heatmap' for now")
    return.plot <- "heatmap"
  }
  if(return.plot == FALSE){
    message("Please set return.plot to either 'volcano', 'heatmap' or 'none'")
    message("Setting return.plot to 'none' for now")
    return.plot <- "none"
  }
  if(missing(label.n.top.genes)) label.n.top.genes <- 30
  if(missing(genes.to.label)) genes.to.label <- "none"
  if(missing(logFC)) logFC <- TRUE
  if(missing(up.genes)) up.genes <- 0
  if(missing(down.genes)) down.genes <- 0
  if(missing(base)){base <- 2}
  if(missing(min.pct)) min.pct <- 0.25
  platypus.version <- "does not matter"
  if(missing(color.p.threshold)){color.p.threshold <- 0.01}
  if(missing(color.log.threshold)){color.log.threshold <- 0.25}
  if(missing(color.by.threshold)){color.by.threshold <- TRUE}
  if(missing(size.top.colorbar)) {size.top.colorbar <- 5}

  if (missing(filter)) filter <- c("MT-", "RPL", "RPS")

  if(!inherits(FindMarkers.out, "data.frame")){

    if(missing(grouping.column)) grouping.column <- "sample_id"

    if(grouping.column %in% names(GEX@meta.data)){
      Seurat::Idents(GEX) <- GEX@meta.data[,c(grouping.column)]
    } else {
      stop("Please provide a valid GEX@meta.data column name as grouping.column")
    }


    cluster_markers <- Seurat::FindMarkers(GEX, min.pct = min.pct, ident.1 = as.character(group1), ident.2 = as.character(group2), base=base)

  } else if (inherits(FindMarkers.out, "data.frame")){
    cluster_markers <- FindMarkers.out
    group1 <- "1"
    group2 <- "2"
  }

  colnames(cluster_markers)[2] <- "avg_logFC"
  SYMBOL <- NULL
  cluster_markers$SYMBOL <- rownames(cluster_markers)

  #adding columns for coloring by threshold
  cluster_markers$signif_pval <- FALSE
  cluster_markers$signif_pval[cluster_markers$p_val_adj < color.p.threshold] <- TRUE
  cluster_markers$signif_logfc <- FALSE
  cluster_markers$signif_logfc[abs(cluster_markers$avg_logFC) > color.log.threshold] <- TRUE
  cluster_markers$color_threshold <- FALSE #adding a column that combines both p and log thresholds (and gate)
  #Just consider pval
  cluster_markers$color_threshold[cluster_markers$signif_pval == TRUE] <- TRUE

  significant_color = cluster_markers$color_threshold
  color_all = all(significant_color,TRUE)

  exclude <- c()
  for (j in filter) {
    exclude <- c(exclude, stringr::str_which(rownames(cluster_markers), j))
  }
  if(length(exclude) > 0){
    cluster_markers <- cluster_markers[-exclude,]
  }
  if (return.plot=="heatmap") {
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
    plot.out <- Seurat::DoHeatmap(GEX, features = heatmap_genes, label = TRUE, size = size.top.colorbar) #Labels turned off as a temporary fix for Seurat related issue. See https://github.com/satijalab/seurat/issues/2399
  } else if (return.plot == "volcano"){

    if (logFC==TRUE) {
      ranks <- order(-cluster_markers$avg_logFC)
      cluster_markers <- cluster_markers[ranks,]
    }
    if (logFC==FALSE) {
      ranks <- order(cluster_markers$p_val_adj)
      cluster_markers <- cluster_markers[ranks,]
    }

    #choose which points to label

    cluster_markers_2 = cluster_markers
    cluster_markers_2 = cluster_markers_2[which(cluster_markers_2$p_val_adj<=0.05),]

    if  (down.genes+up.genes != 0) {

      genes_to_label_asc = data.frame()
      genes_to_label_desc = data.frame()

      if (down.genes>0) {
        cluster_markers_asc <- cluster_markers_2[order(cluster_markers_2$avg_logFC),] #From - to + Expresssion
        genes_to_label_asc = cluster_markers_asc[1:down.genes,]
      }

      if (up.genes>0) {
        cluster_markers_desc <- cluster_markers_2[order(cluster_markers_2$avg_logFC, decreasing = TRUE),]   #From + to - Expression descending
        genes_to_label_desc <- cluster_markers_desc[1:up.genes,]
      }

      cluster_markers_rel <- rbind(genes_to_label_asc,genes_to_label_desc )

    } else if  (down.genes+up.genes == 0){

      cluster_markers_rel <- cluster_markers[1:label.n.top.genes,]
    }

    if(genes.to.label[1] != "none"){
      extra_genes <- subset(cluster_markers, SYMBOL %in% genes.to.label)
      if(nrow(extra_genes) > 0){
        cluster_markers_rel <- rbind(cluster_markers_rel, extra_genes)
      }
    }

    if(color.by.threshold == FALSE){
      plot.out <- ggplot2::ggplot(cluster_markers, ggplot2::aes(x = avg_logFC, y = -log10(p_val_adj), col = avg_logFC)) + ggplot2::geom_point(show.legend = F, size = 3, alpha = 0.7) + cowplot::theme_cowplot() + ggplot2::theme(panel.background = ggplot2::element_blank(),axis.text = ggplot2::element_text(size = 30), axis.line = ggplot2::element_line(size = 2), axis.ticks = ggplot2::element_line(size = 2), axis.ticks.length = ggplot2::unit(0.3, "cm"), text = ggplot2::element_text(size=30), legend.position = "none") + ggplot2::labs(title = paste0("DEGs ", group1, " vs. ", group2), x = "log2(FC)", y = "-log10(adj p)") + ggrepel::geom_text_repel(data = cluster_markers_rel, ggplot2::aes(x = avg_logFC, y = -log10(p_val_adj), label = SYMBOL), inherit.aes = F, size = 6, segment.alpha = 1, max.overlaps = 50) + ggplot2::scale_colour_viridis_c(option = "B")

    } else if(color.by.threshold == TRUE && color_all == FALSE){

      plot.out <- ggplot2::ggplot(cluster_markers, ggplot2::aes(x = avg_logFC, y = -log10(p_val_adj), col = color_threshold)) + ggplot2::geom_point(show.legend = F, size = 3, alpha = 0.7) + cowplot::theme_cowplot() + ggplot2::theme(panel.background = ggplot2::element_blank(),axis.text = ggplot2::element_text(size = 30), axis.line = ggplot2::element_line(size = 2), axis.ticks = ggplot2::element_line(size = 2), axis.ticks.length = ggplot2::unit(0.3, "cm"), text = ggplot2::element_text(size=30), legend.position = "none") + ggplot2::labs(title = paste0("DEGs ", group1, " vs. ", group2), x = "log2(FC)", y = "-log10(adj p)") + ggrepel::geom_text_repel(data = cluster_markers_rel, ggplot2::aes(x = avg_logFC, y = -log10(p_val_adj), label = SYMBOL), inherit.aes = F, size = 6, segment.alpha = 1, max.overlaps = 50) + ggplot2::scale_colour_manual(values = c("black","darkred"))

    } else if(color.by.threshold == TRUE && color_all == TRUE){

      plot.out <- ggplot2::ggplot(cluster_markers, ggplot2::aes(x = avg_logFC, y = -log10(p_val_adj), col = color_threshold)) + ggplot2::geom_point(show.legend = F, size = 3, alpha = 0.7) + cowplot::theme_cowplot() + ggplot2::theme(panel.background = ggplot2::element_blank(),axis.text = ggplot2::element_text(size = 30), axis.line = ggplot2::element_line(size = 2), axis.ticks = ggplot2::element_line(size = 2), axis.ticks.length = ggplot2::unit(0.3, "cm"), text = ggplot2::element_text(size=30), legend.position = "none") + ggplot2::labs(title = paste0("DEGs ", group1, " vs. ", group2), x = "log2(FC)", y = "-log10(adj p)") + ggrepel::geom_text_repel(data = cluster_markers_rel, ggplot2::aes(x = avg_logFC, y = -log10(p_val_adj), label = SYMBOL), inherit.aes = F, size = 6, segment.alpha = 1, max.overlaps = 50) + ggplot2::scale_colour_manual(values = c("darkred"))

    }
  }
  if (return.plot=="none") plot.out <- NULL
  return(list(cluster_markers, plot.out))
}
