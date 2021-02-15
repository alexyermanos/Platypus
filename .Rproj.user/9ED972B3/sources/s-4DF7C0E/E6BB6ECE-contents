#' Integrates VDJ and gene expression libraries by providing cluster membership seq_per_vdj object and the index of the cell in the Seurat RNA-seq object.
#' @param automate_GEX.output A single seurat object from automate_GEX function. This will likely be supplied as automate_GEX.output[[1]] for an individual set of integrated repertoires.
#' @param clonotype.ids  Output from either VDJ_analyze or VDJ_clonotype functions. This list should correspond to a single GEX.list object, in which each list element in clonotype.list is found in the GEX.object. Furthermore, these repertoires should be found in the automate_GEX library.
#' @return Returns a stacked barplot that visualizes the seurat cluster membership for different cell phenotypes.
#' @export
#' @examples
#' \dontrun{
#' GEX_phenotype_per_clone_plot <- GEX_phenotype_per_clone(GEX.object = automate.gex.output[[1]], clonotype.ids= c(1,2,3,4,5))
#'}
#'
GEX_phenotype_per_clone <- function(automate_GEX.output,clonotype.ids){
  require(ggplot2)
  require(reshape2)
  dca<-dcast(seurat.object@meta.data,clonotype_id~cell.state,value.var="cell.state",length)
  dca$Sum <- rowSums(dca[-1])
  dca[2:(ncol(dca)-1)]<- dca[2:(ncol(dca)-1)]/dca$Sum
  dca<-melt(dca[-ncol(dca)], id.vars = "clonotype_id")
  dca<-dca[which(dca$clonotype_id %in% clonotype.ids),]
  stacked.ggplot<-ggplot(data=dca, aes(x=clonotype_id, y=value, fill=variable)) +
    geom_bar(stat="identity")+
    ylab("Cell Counts")+
    xlab("Clonotypes")+
    labs(fill = "Cell State")
  return(stacked.ggplot)
}
