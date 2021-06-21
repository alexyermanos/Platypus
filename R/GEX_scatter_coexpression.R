#' Clonal frequency plot displaying clonal expansion for either T and B cells with Platypus v3 input. 
#' @param VDJ.matrix GEX seurat object generated with VDJ_GEX_matrix
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
for_dot <- FetchData(GEX.matrix, vars = c(gene.1,gene.2))


min_x <- round(range(for_dot[,1])[1] - 0.15,1)
min_y <- round(range(for_dot[,2])[1] - 0.15,1)
max_x <- round(range(for_dot[,1])[2] + 0.15,1)
max_y <- round(range(for_dot[,2])[2] + 0.15,1)

hist_x <- ggplot(data = for_dot,aes(x = for_dot[,1])) +geom_density(adjust = 0.4, alpha = 0.2, col = color.theme, fill = color.theme, size = 2) +theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank(), axis.ticks.x = element_blank() , axis.text.x = element_blank(),panel.background = element_blank(),panel.border = element_blank(),text = element_text(size = 30), plot.margin = unit(c(3,-5.5,4,6.5), "mm"),axis.line.y.left = element_line(size = 1.2)) + xlim(c(min_x, max_x)) 

hist_y <- ggplot(data = for_dot,aes(y = for_dot[,2])) +geom_density(adjust = 0.4, alpha = 0.2, col = color.theme, fill = color.theme, size = 2)+theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x = element_text(color = "white"),panel.background = element_blank(),text = element_text(size = 30), plot.margin = unit(c(3,-5.5,4.1,3), "mm"),axis.line.x.bottom = element_line(size = 1.2)) +labs(x = "") + ylim(c(min_y, max_y))

empty <- ggplot() +geom_point(aes(1,1), colour="white")+ theme(axis.ticks=element_blank(), panel.background=element_rect(fill = "transparent"),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.title.y=element_blank())

scatter <- ggplot(for_dot, aes(x = for_dot[,1], y = for_dot[,2], color = (for_dot[,1] + for_dot[,2]))) + geom_point(size = 3) + theme(plot.margin = unit(c(3,-5.5,4,3), "mm"),panel.background = element_blank(),text = element_text(size = 30), legend.position = "none", axis.line = element_line(size = 1.2)) + xlim(c(min_x, max_x)) + ylim(c(min_y, max_y))+ labs(x = gene.1, y = gene.2) + scale_colour_gradient2(low = "grey40", mid = "grey70", high = color.theme)


return(grid.arrange(hist_x, empty , scatter, hist_y, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4)))

}

