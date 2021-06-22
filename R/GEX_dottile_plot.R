#' Outputs a dotplot for gene expression, where the color of each dot is scaled by the gene expression level and the size is scaled by the \% of cells positive for the gene
#' @param GEX.matrix GEX seurat object generated with VDJ_GEX_matrix
#' @param genes Character vector. Genes of those in rownames(GEX.matrix) to plot. Can be any number, but more then 30 is discuraged because of cluttering
#' @param group.by Character. Name of a column in GEX.matrix@meta.data to split the plot by. If set to "none", a plot with a single column will be produced.
#' @threshold.to.plot Integer 1-100. % of cells which must be expressing the feature to plot a point. If below, the field will be left empty
#' @param platypus.version This is coded for "v3" only, but in practice any Seurat Object can be fed in
#' @return Returns a ggplot object
#' @export
#' @examples
#' \dontrun{
#' To return a plot detailing the expression of common genes by seurat cluster
#'GEX_dottile_plot(GEX.matrix = VDJ.GEX.matrix.output[[2]], genes = c("CD19", "EBF1","CD3E","CD8A","MKI67"), group.by = "seurat_clusters", threshold.to.plot = 5)
#'}
GEX_dottile_plot <- function(GEX.matrix,
                             genes,
                             group.by,
                             threshold.to.plot,
                             platypus.version){
  platypus.version <- "v3"

  if(missing(threshold.to.plot)) threshold.to.plot <- 5

  if(missing(genes)){
    if(GEX.matrix$celltype[1] == "T cell"){
      genes <- c("CD4","CD8A","CD28","ICOS","CD40LG","PDCD1","LAG3","S1PR1","PTPRC","SELL","TBX21","GATA3","SPI1","IRF4","BCL6","STAT4","STAT6","FOXP3", "MKI67","TCF7","KLRG1","TOX", "GZMB","IFNG","TGFB1","IL10","IL17A","CSF2","IL2RA","IL4RA","IL12RB1","CXCR3","CXCR5","CCR5","CCR8","SELL","ITGA4","ITGB2")
    } else {
      genes <- c("CD19", "CD74","SDC1", "EBF1","PTPRC","CD93","CD38","CD24A","CD34","CD1D1","CR2","MS4A1","CXCR5","SELL","CD40","CD83","H2-AB1","H2-EB1","CD27","POU2AF1","NT5E","FAS","PDCD1LG2","PRDM1","ITGAM","IL10","IL12A","HAVCR2")
    }
  }

  require(tidyr)
  if(missing(GEX.matrix)) stop("Please provide a seurat object as input to GEX.matrix")

  to_del <- c()
  for(i in 1:length(genes)){
    if(!genes[i] %in% rownames(GEX.matrix)){
      print(paste0(genes[i], " not found in seurat object. This gene is skipped"))
      to_del <- c(to_del, i)
    }
  }
  if(length(to_del) > 0){
    genes <- genes[-to_del]
  }

  if(missing(group.by)) group.by <- "none"
  if(!group.by %in% names(GEX.matrix@meta.data)){
    print("group.by column not found. Returning plot with single column")
    group.by <- "singlegroup"
    GEX.matrix$singlegroup <- "group"
  }


  to_plot_f <- FetchData(GEX.matrix, vars = c(group.by, genes))

  #pivot
  to_plot_f <- pivot_longer(to_plot_f, cols = c(2:ncol(to_plot_f)))
  names(to_plot_f)[1] <- "group"
  #summarize and group
  #for now by cluster
  to_plot_sum_f <- to_plot_f %>% group_by(group,name) %>% summarise(mean_scaled_expression = mean(value[value > 0]), perc_expressing_cells = (length(value[value > 0])/n()*100))

  to_plot_sum_f$name <- ordered(as.factor(to_plot_sum_f$name), levels = rev(genes))

  to_plot_sum_f$perc_expressing_cells[to_plot_sum_f$perc_expressing_cells < threshold.to.plot] <- NA

  plot_out <- ggplot(to_plot_sum_f, aes(x = group, y = name, col = mean_scaled_expression, size = perc_expressing_cells)) + geom_point(show.legend = T)  + theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=3),axis.text = element_text(size = 30), axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1), axis.line = element_blank(), axis.ticks = element_line(size = 2), axis.ticks.length = unit(0.3, "cm"), text = element_text(size=30), legend.key = element_rect(fill = "white", color = "white"), legend.position = "right") + labs(title = paste0("Expression by ", group.by), x = "", y = "", color = "Scaled expression", size = "% of expressing cells")  + scale_color_viridis_c(option = "B", end = 0.9) +scale_size_binned(range = c(1,9.5))

    return(plot_out)
}

