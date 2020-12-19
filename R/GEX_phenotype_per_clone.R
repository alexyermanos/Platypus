GEX_phenotype_per_clone <- function(seurat.object,clonotype.ids){
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