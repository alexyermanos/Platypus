#' Visualization of marker expression in a data set or of predefined genes (B cells, CD4 T cells and CD8 T cells).
#' @param GEX GEX output of the VDJ_GEX_matrix function (VDJ_GEX_matrix[[2]])).
#' @param gene_set Character vector containing the markers of interest given by the user.
#' @param predefined_genes Character vector to chose between B_cell, CD4_T_cell, and CD8_T_cell.
#' @return Return a list. Element[[1]] is the feature plot of markers of interest or predefined genes. Element[[2]] is the dottile plot of markers of interest or predefined genes. Element[[3]] is the violin plot of markers of interest or predefined genes.
#' @export
#' @examples
#' \dontrun{
#' # Pre-defined gene set for CD4 T cells
#' GEX_gene_visualization(GEX = VGM$GEX, predefined_genes = "CD4_T_cell")
#' 
#' # Pre-defined gene set for CD8 T cells
#' GEX_gene_visualization(GEX = VGM$GEX, predefined_genes = "CD8_T_cell")
#' 
#' # Pre-defined gene set for B cells
#' GEX_gene_visualization(GEX = VGM$GEX, predefined_genes = "B_cell")
#' 
#' # Gene set defined by user
#' GEX_gene_visualization(GEX = VGM$GEX, gene_set=c("CD8A","CD3E","SELL","FAS","ID3","SDC1"))
#'}

GEX_gene_visualization<-function(GEX, gene_set, predefined_genes=c("B_cell","CD4_T_cell","CD8_T_cell"), group.by){
  if(missing(GEX))stop("Please provide GEX input for this function")
  if(missing(group.by)){
    group.by = "sample_id"
  }
  platypus.version <- "It doesn't matter"
  
  if(missing(gene_set)){
    if(missing(predefined_genes))stop("Please provide a gene set or use predefined_genes (B_cell, CD4_T_cell or CD8_T_cell) input for this function")
    if(predefined_genes == "B_cell"){
      genes<-c("SLPI","XBP1","SLAMF7","SDC1","PRDM1","H2-AA","H2-EB1","H2-EB2","H2-OA","H2-OB","IL10RA","IFNAR1","IFNAR2","IL21R","CXCR3","CXCR4","CXCR5","CCR6","CCR7","MS4A1","IRF8","DUSP1","CD19","EBF1","PAX5","TCF3","SELL","FAS","NR4A1","CR2","FOXP1","CD38","CD69","PTPRC","CD79A","AICDA")
    } else if (predefined_genes=="CD4_T_cell"){
      genes<-c("ID3","CCR7","IL7R","LEF1","SLAMF6","BCL2","TCF7","CXCR6","TBX21","ID2","GZMB","IFNG","LY6C2","LAG3","CXCR5","BCL6","ASCL2","PDCD1","ICOS","IL21","IL4","FOXP3","IL2RA","IL10","CD81","CD74","KLRG1")
    } else if (predefined_genes=="CD8_T_cell"){
      genes<-c("SELL","CD69","TIGIT","IFNG","GZMB","GZMA","TBX21","NFATC1","SLAMF6","EOMES","ID2","STAT3","STAT4","FOXO1","CD3E", "CD8A","CD44")
    }
  } else{
    genes<-gene_set
  }
  #Feature plot
  feature_plot<-FeaturePlot(GEX, features = genes)
  #Dottile function
  Dottile_plot<-GEX_dottile_plot(GEX, genes = genes, group.by = group.by)
  #Violin plot for markers
  violin_plot<-VlnPlot(GEX, features = genes)
  results<-list()
  results$feature_plot<-feature_plot
  results$dottile_plot<-Dottile_plot
  results$violin_plot<-violin_plot
  return(results)
}
