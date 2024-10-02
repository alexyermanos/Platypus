#'Coexpression of selected genes
#'
#'@description Returns eiter a plot or numeric data of coexpression levels of selected genes.Coexpression \% is calculated as the quotient of double positive cells (counts \> 0) and the sum of total cells positive for either genes.
#' @param GEX GEX seurat object generated with VDJ_GEX_matrix (VDJ_GEX_matrix.output\[\[2\]\])
#' @param genes Character vector. At least 2 genes present in rownames(GEX). Use "all" to include all genes. The number of comparisons to make is the length(genes)! (factorial). More than 100 genes are not recommended.
#' @param subsample.n Interger. Number of cells to subsample. If set to 100, 100 cells will be randomly sampled for the calculation
#' @param plot.dotmap Boolean. Whether to return a plot
#' @return Returns a dataframe if pot.dotmap == FALSE or a ggplot if plot.dotmap == TRUE detailing the coexpression levels of selected genes within the given cell population
#' @export
#' @examples
#' GEX_coexpression_coefficient(GEX = Platypus::small_vgm[[2]]
#', genes = c("CD19", "CD83"), subsample.n = "none", plot.dotmap = FALSE)


GEX_coexpression_coefficient <- function(GEX,
                                         genes,
                                         subsample.n,
                                         plot.dotmap){
  gene.1 <- NULL
  gene.2 <- NULL
  coex.coef <- NULL
  perc.single.positive <- NULL

  if(missing(plot.dotmap)) plot.dotmap <- TRUE

  coex_coef <- function(combs,cmat){
    #Nr of douple positives / nr of single positives of both markers
    return(sum(((cmat[,combs[1]] > 0) + (cmat[,combs[2]] > 0)) == 2)  / (sum(cmat[,combs[1]] > 0) + sum(cmat[,combs[2]] > 0)) * 100)
  }
  coex_n <- function(combs,cmat){
    #Nr of nr of single positives of both markers
    return((sum(c(cmat[,combs[1]] > 0),(cmat[,combs[2]] > 0))) / 2*nrow(cmat) * 100)
  }

  if(missing(genes)){genes <- "all"}
  if(missing(subsample.n)) subsample.n <- "none"

  if(genes[1] != "all" & length(genes > 1)){
    cmat <- SeuratObject::FetchData(GEX, vars = genes, slot = "counts")
  } else {
    #get all counts
    cmat <- GEX@assays$RNA@counts
  }

  if(subsample.n[1] != "none"){
    cmat <- cmat[sample(1:nrow(cmat),subsample.n),]
  }

  combs <- utils::combn(c(1:ncol(cmat)), m = 2, simplify = FALSE)

  message(paste0("Calculating coexpression for ", ncol(cmat), " genes with ", nrow(cmat), " cells"))

  out_coef <- do.call("c", lapply(combs, coex_coef ,cmat))
  out_n_dp <- do.call("c", lapply(combs, coex_coef ,cmat))

  col.1 <- do.call("c",lapply(combs, function(x) return(x[1])))
  col.2 <- do.call("c",lapply(combs, function(x) return(x[2])))

  out_t <- data.frame("col.1" = col.1,
                      "col.2" = col.2,
                      "gene.1" = colnames(cmat)[col.1],
                      "gene.2" = colnames(cmat)[col.2],
                      "coex.coef" = out_coef,
                      "perc.single.positive" = out_n_dp)

  out_t$gene.1 <- as.factor(out_t$gene.1)
  out_t$gene.2 <- as.factor(out_t$gene.2)
  out_t <- out_t[order(out_t$gene.1),]
  out_t <- out_t[order(out_t$gene.2),]
  out_t$gene.1 <- ordered(out_t$gene.1, levels = genes)
  out_t$gene.2 <- ordered(out_t$gene.2, levels = genes)

  if(plot.dotmap == TRUE){

    plot_out <- ggplot2::ggplot(out_t, ggplot2::aes(x = gene.2, y = gene.1, col = coex.coef, size = perc.single.positive)) + ggplot2::geom_point(show.legend = TRUE) + cowplot::theme_cowplot()  + ggplot2::labs(title = "", x = "", y = "", color = "Coexpression coefficient", size = "% of single positives") + ggplot2::scale_color_viridis_c(option = "B", end = 0.9) + ggplot2::scale_x_discrete(limits = rev)

    #+ ggplot2::scale_size_binned(breaks = c(5,25,50,75), labels = c("<5","5", "25",">50"),range = c(2,7))
    #panel.background = ggplot2::element_blank(),axis.text = ggplot2::element_text(size = 30), axis.text.x = ggplot2::element_text(angle = 60, vjust = 0.95, hjust=1), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_line(size = 2), axis.ticks.length = ggplot2::unit(0.3, "cm"), text = ggplot2::element_text(size=30), legend.key = ggplot2::element_rect(fill = "white", color = "white"),

    return(plot_out)

  } else {
    return(out_t)
  }

}
