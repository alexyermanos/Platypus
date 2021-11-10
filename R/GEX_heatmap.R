#' Produces a heatmap containing gene expression information at the clonotype level. The rows correspond to different genes that can either be determined by pre-made sets of B or T cell markers, or can be customized by the user. The columns correspond to individual cells and the colors correspond to the different clonotype families.
#' @param GEX A single seurat object from clonotype_GEX function corresponding to all of the samples in a single VDJ_analyze object. This will likely be supplied as clonotype_GEX.output[[i]] if there were multiple, distinct transcriptomes.
#' @param b.or.t Logical indicating if B or T cell gene panel should be used.
#' @param sample.index Corresponds to which repertoire should be used in the case that the length of clonotype.list has a length greater than 1. The transcriptional profiles from only one repertoire can be plotted at a time.
#' @param clone.rank.threshold A numeric that specifies the threshold clonal rank that specifies which clonotypes to extract transcriptome information from. For example, if 10 is supplied then the gene expression for the top ten clones included on the heatmap, separated by clonotype.
#' @param custom.array Corresponds to which repertoire should be used in the case that the length of clonotype.list has a length greater than 1. The transcriptional profiles from only one repertoire can be plotted at a time.
#' @param slot Seurat data slot from which to plot values. Can be "raw.data", "data" or "scale.data"
#' @return Returns a heatmap via Seurat::DoHeatmap of gene expression per clonotype
#' @seealso VDJ_extract_sequences
#' @export
#' @examples
#' #prep the small_vgm sample dataset
#' small_vgm <- Platypus::small_vgm
#' small_vgm[[2]]$clone_rank <- c(1:nrow(small_vgm[[2]]@meta.data))
#' GEX_heatmap(GEX = small_vgm[[2]],b.or.t = "custom"
#' ,clone.rank.threshold = 1,sample.index = "s1"
#' ,custom.array = c("CD24A","CD83"), slot = "data")
#'

GEX_heatmap <- function(GEX,b.or.t,sample.index,clone.rank.threshold,custom.array, slot){
  GEX.object <- GEX
  if(missing(b.or.t)) b.or.t <- "t"
  if(missing(slot)) slot <- "scale.data"
  if(missing(custom.array)) custom.array <- c("")
  holding_sample_id <- which(GEX.object$sample_id==sample.index)
  holding_clone_rank_index <- which(GEX.object$clone_rank<=clone.rank.threshold)
  holding_cells <- intersect(holding_sample_id,holding_clone_rank_index)
  if(b.or.t=="b"){
    GEX.heatmap <- Seurat::DoHeatmap(GEX.object,features = c("CD74","CD79A","FCER2","TMSB10","BTG1","BACH2","MEF2C","HVCN1","SARAF","CXCR4","FCRL1","CD72","NCF1","AIM2","CRIP1","CD82","ITGB1","CD24","PTPRC","CD19","CD27","CD38","SDC1","CD22","FAS","TNFRSF13B","IL4R","XBP1"),cells = holding_cells,label = F,group.by = "clone_rank",lines.width = 1, slot = slot)
  }
  else if(b.or.t=="t"){
    GEX.heatmap <- Seurat::DoHeatmap(GEX.object,features = c("CD3E","CD8A","CD4","IL7R","CCL5","GZMA","GZMK","GZMB","PRF1","NKG7","SELL","TNFRSF13B","TBX1","PDCD1","ITGAE","LAG3","CD44","ICOS","CX3CR1","MKI67","TCF7","CST7","GNLY","CXCR5","EOMES"),cells = holding_cells,label = F,group.by = "clone_rank",lines.width = 1, slot = slot)
  }
  else if(b.or.t=="custom"){
    GEX.heatmap <- Seurat::DoHeatmap(GEX.object,features = custom.array,cells = holding_cells,label = F,group.by = "clone_rank",lines.width = 1, slot = slot)
  }
  return(GEX.heatmap)
}
