#' Clonal frequency plot displaying clonal expansion for either T and B cells with Platypus v3 input.
#' @param GEX.matrix GEX seurat object generated with VDJ_GEX_matrix
#' @param gene.1 Character. Name of a gene in rownames(VDJ.matrix)
#' @param gene.2 Character. Name of a gene in rownames(VDJ.matrix)
#' @param color.theme Character. A color to use for the composite plot
#' @return Returns a gridplot showing coexpression of gene.1 and gene.2
#' @export
#' @examples
#' \dontrun{
#'
#'gene1 <- "BCL6"
#'gene2 <- "CD40"
#'ggsave(filename = paste0("Coex_", gene1,"_", gene2,"_scat.png"), plot = GEX_scatter_coexpression(VDJ_comb[[2]], gene1,gene2), dpi = 300, width = 12, height = 12)
#'
#'}

GEX_scatter_coexpression <- function(GEX.matrix,
                                     gene.1,
                                     gene.2,
                                     color.theme){
require(gridExtra)

if(missing(color.theme)) color.theme <- "darkorchid4"

#get data
for_dot <- SeuratObject::FetchData(GEX.matrix, vars = c(gene.1,gene.2))


min_x <- round(range(for_dot[,1])[1] - 0.15,1)
min_y <- round(range(for_dot[,2])[1] - 0.15,1)
max_x <- round(range(for_dot[,1])[2] + 0.15,1)
max_y <- round(range(for_dot[,2])[2] + 0.15,1)

hist_x <- ggplot2::ggplot(data = for_dot,ggplot2::aes(x = for_dot[,1])) + ggplot2::geom_density(adjust = 0.4, alpha = 0.2, col = color.theme, fill = color.theme, size = 2) + ggplot2::theme(legend.position = "none",axis.title.x = ggplot2::element_blank(),axis.title.y = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank() , axis.text.x = ggplot2::element_blank(),panel.background = ggplot2::element_blank(),panel.border = ggplot2::element_blank(),text = ggplot2::element_text(size = 30), plot.margin = ggplot2::unit(c(3,-5.5,4,6.5), "mm"),axis.line.y.left = ggplot2::element_line(size = 1.2)) + ggplot2::xlim(c(min_x, max_x))

hist_y <- ggplot2::ggplot(data = for_dot, ggplot2::aes(y = for_dot[,2])) + ggplot2::geom_density(adjust = 0.4, alpha = 0.2, col = color.theme, fill = color.theme, size = 2)+ ggplot2::theme(legend.position = "none",axis.title.y = ggplot2::element_blank(),axis.text.y = ggplot2::element_blank() , axis.ticks.y = ggplot2::element_blank(),axis.title.x = ggplot2::element_text(color = "white"),panel.background = ggplot2::element_blank(),text = ggplot2::element_text(size = 30), plot.margin = ggplot2::unit(c(3,-5.5,4.1,3), "mm"),axis.line.x.bottom = ggplot2::element_line(size = 1.2)) +ggplot2::labs(x = "") + ggplot2::ylim(c(min_y, max_y))

empty <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(1,1), colour="white")+ ggplot2::theme(axis.ticks= ggplot2::element_blank(), panel.background= ggplot2::element_rect(fill = "transparent"),axis.text.x= ggplot2::element_blank(), axis.text.y= ggplot2::element_blank(),axis.title.x= ggplot2::element_blank(), axis.title.y= ggplot2::element_blank())

scatter <- ggplot2::ggplot(for_dot, ggplot2::aes(x = for_dot[,1], y = for_dot[,2], color = (for_dot[,1] + for_dot[,2]))) + ggplot2::geom_point(size = 3) + ggplot2::theme(plot.margin = ggplot2::unit(c(3,-5.5,4,3), "mm"),panel.background = ggplot2::element_blank(),text = ggplot2::element_text(size = 30), legend.position = "none", axis.line = ggplot2::element_line(size = 1.2)) + ggplot2::xlim(c(min_x, max_x)) + ggplot2::ylim(c(min_y, max_y))+ ggplot2::labs(x = gene.1, y = gene.2) + ggplot2::scale_colour_gradient2(low = "grey40", mid = "grey70", high = color.theme)


return(gridExtra::grid.arrange(hist_x, empty , scatter, hist_y, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4)))
}

