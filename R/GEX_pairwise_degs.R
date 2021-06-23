#' Produces and saves a list of volcano plots with each showing differentially expressed genes between pairs groups. If e.g. seurat_clusters used as group.by, a plot will be generated for every pairwise comparison of clusters. For large numbers of this may take longer to run. Only available for platypus v3
#' @param GEX Output Seurat object of the VDJ_GEX_matrix function (VDJ_GEX_matrix.output[[2]])
#' @param group.by Character. Defaults to "seurat_clusters" Column name of GEX@meta.data to use for pairwise comparisons. More than 20 groups are discuraged.
#' @param min.pct Numeric. Defaults to 0.25 passed to Seurat::FindMarkers
#' @param RP.MT.filter Boolean. Defaults to True. If True, mitochondrial and ribosomal genes are filtered out from the output of Seurat::FindMarkers
#' @param label.n.top.genes Integer. Defaults to 50. Defines how many genes are labelled via geom_text_repel. Genes are ordered by adjusted p value and the first label.n.genes are labelled
#' @param genes.to.label Character vector. Defaults to "none". Vector of gene names to plot indipendently of their p value. Can be used in combination with label.n.genes.
#' @param save.plot Boolean. Defaults to True. Whether to save plots as appropriately named .png files
#' @return A nested list with out[[i]][[1]] being plots and out[[i]][[2]] being DEG dataframes.
#' @export
#' @examples
#' \dontrun{
#' GEX_pairwise_DEGs <- function(GEX = VDJ_GEX_matrix.output[[2]],group.by = "seurat_clusters"
#' ,min.pct = 0.25,RP.MT.filter = T,label.n.top.genes = 50,genes.to.label = c("DIABLO","ELMO1")
#' ,save.plot = F)
#'}
GEX_pairwise_DEGs <- function(GEX,
                              group.by,
                              min.pct,
                              RP.MT.filter,
                              label.n.top.genes,
                              genes.to.label,
                              save.plot){

  gene <- NULL
  avg_log2FC <- NULL
  p_val_adj <- NULL
  perc_expressing_cells <- NULL

  require(ggrepel)

  if(missing(label.n.top.genes)) label.n.top.genes <- 50
  if(missing(genes.to.label)) genes.to.label <- "none"
  if(missing(group.by)) group.by <- "sample_id"
  if(missing(min.pct)) min.pct <- 0.25
  if(missing(RP.MT.filter)) RP.MT.filter <- T
  if(missing(save.plot)) save.plot <- F

  if(!group.by %in% names(GEX@meta.data)){stop("Please enter valid metadata column name")}

  to_group <- unique(as.character(GEX@meta.data[,group.by]))
  SeuratObject::Idents(GEX) <- as.character(GEX@meta.data[,group.by])

  if(length(to_group) == 1){stop("Grouping column has to contain at least two unique entries")}

  if(length(to_group) > 2){
    combs <- as.data.frame(t(utils::combn(to_group, m = 2,simplify = TRUE)))#get combinations to test

    combs[,1] <- ordered(as.factor(combs[,1]), levels = rev(to_group))
    combs[,2] <- ordered(as.factor(combs[,2]), levels = to_group)

  } else {
    combs <- data.frame(to_group[1], to_group[2])
  }

  degs.list <- list()
  plot.list <- list()
  for(i in 1:nrow(combs)){
    print(paste0("Calculating pairwise DEGs ", i, " of ", nrow(combs)))
    print(combs[i,1])
    print(combs[i,2])
    degs <- Seurat::FindMarkers(GEX, ident.1 = combs[i,1], ident.2 = combs[i,2],min.pct = min.pct)

    degs$gene <- rownames(degs)

    if(RP.MT.filter == T){
      degs <- subset(degs, stringr::str_detect(degs$gene, "(^RPL)|(^MRPL)|(^MT-)|(^RPS)") == F)
    }

    #choose which points to label
      degs <- degs[order(degs$p_val_adj),]
      degs_rel <- degs[1:label.n.top.genes,]

      if(genes.to.label[1] != "none"){
        extra_genes <- subset(degs, gene %in% genes.to.label)
        if(nrow(extra_genes) > 0){
          degs_rel <- rbind(degs_rel, extra_genes)
        }
      }

    plot.out <- ggplot2::ggplot(degs, ggplot2::aes(x = avg_log2FC, y = -log10(p_val_adj), col = avg_log2FC)) + ggplot2::geom_point(show.legend = F, size = 3, alpha = 0.7) + ggplot2::theme(legend.position = "none",panel.background = ggplot2::element_blank(),axis.text = ggplot2::element_text(size = 30), axis.line = ggplot2::element_line(size = 2), axis.ticks = ggplot2::element_line(size = 2), axis.ticks.length = ggplot2::unit(0.3, "cm"), text = ggplot2::element_text(size=30)) + ggplot2::labs(title = paste0("DEGs ", combs[i,1], " vs. ", combs[i,2]), x = "log2(FC)", y = "-log10(adj p)")+ ggrepel::geom_text_repel(data = degs_rel, ggplot2::aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), inherit.aes = F, size = 6, segment.alpha = 1, max.overlaps = 50) + ggplot2::scale_colour_viridis_c(option = "B")

    if(save.plot == T){
      ggplot2::ggsave(plot.out, filename = paste0("DEGs_", combs[i,1], "_vs_", combs[i,2],".png"), dpi = 400, width = 10, height = 10)
    }
    degs.list[[i]] <- degs
    plot.list[[i]] <- plot.out
  }
  return(list(plot.list, degs.list))
}

