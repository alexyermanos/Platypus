#' Clonal frequency plot displaying clonal expansion for either T and B cells with Platypus v3 input. 
#' @param GEX.matrix GEX seurat object generated with VDJ_GEX_matrix
#' @param genes Character vector. At least 2 genes present in rownames(GEX.matrix). Use "all" to include all genes. The number of comparisons to make is the length(genes)! (factorial). More than 100 genes are not recommended. 
#' @param subsample.n Interger. Number of cells to subsample. If set to 100, 100 cells will be randomly sampled for the calculation
#' @param plot.dotmap Boolean. Whether to return a plot
#' @return Returns a dataframe if pot.dotmap == F or a ggplot if plot.dotmap == T
#' @export
#' @examples
#' \dontrun{
#' 
#' To return a dataframe with coefficients
#'coef_out <- GEX_coexpression_coefficient(GEX.matrix = VDJ_comb[[2]], genes = c("CD19", "EBF1","SDC1"), subsample.n = "none", plot.dotmap = F)
#'
#'To return a dotplot detailing coexpression and overall expression
#'plot_out <- GEX_coexpression_coefficient(GEX.matrix = VDJ_comb[[2]], genes = c("CD19", "EBF1","SDC1"), subsample.n = "none", plot.dotmap = T)
#'}

GEX_coexpression_coefficient <- function(GEX.matrix, 
                                         genes, 
                                         subsample.n,
                                         plot.dotmap){
  if(missing(plot.dotmap)) plot.dotmap <- T
  
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
    cmat <- FetchData(GEX.matrix, vars = genes, slot = "counts")
  } else {
    #get all counts
    cmat <- GEX.matrix@assays$RNA@counts
  }
  
  if(subsample.n[1] != "none"){
    cmat <- cmat[sample(1:nrow(cmat),subsample.n),]
  }
  
  combs <- combn(c(1:ncol(cmat)), m = 2, simplify = F)
  
  print(paste0("Calculating coexpression for ", ncol(cmat), " genes with ", nrow(cmat), " cells"))
  
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
  
  if(plot.dotmap == T){
    
    plot_out <- ggplot(out_t, aes(x = gene.2, y = gene.1, col = coex.coef, size = perc.single.positive)) + geom_point(show.legend = T)  + theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=3),axis.text = element_text(size = 30), axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1), axis.line = element_blank(), axis.ticks = element_line(size = 2), axis.ticks.length = unit(0.3, "cm"), text = element_text(size=30), legend.key = element_rect(fill = "white", color = "white"), legend.position = "right") + labs(title = "", x = "", y = "", color = "Coexpression coefficient", size = "% of single positives") + scale_color_viridis_c(option = "B", end = 0.9) +scale_size_binned(breaks = c(5,25,50,75), labels = c("<5","5", "25",">50"),range = c(2,7)) + scale_y_discrete(limits = rev)
    plot_out
    
    return(plot_out)
    
  } else {
    return(out_t)
  }

}
