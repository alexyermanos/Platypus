#' Integrates VDJ and gene expression libraries by providing cluster membership seq_per_vdj object and the index of the cell in the Seurat RNA-seq object.
#' @param seurat.object A single seurat object from automate_GEX function after labeling cell phenotypes using the GEX_phenotype function.
#' @param clonotype.ids  Output from either VDJ_analyze or VDJ_clonotype functions. This list should correspond to a single GEX.list object, in which each list element in clonotype.list is found in the GEX.object. Furthermore, these repertoires should be found in the automate_GEX library.
#' @return Returns a stacked barplot that visualizes the seurat cluster membership for different cell phenotypes.
#' @export
#' @examples
#' \dontrun{
#' GEX_phenotype_per_clone_plot <- GEX_phenotype_per_clone(seurat.object = automate.gex.output[[1]], clonotype.ids= c(1,2,3,4,5))
#'}
#'
GEX_phenotype_per_clone <- function(seurat.object,clonotype.ids){
  require(ggplot2)
  require(reshape2)
  # dca<-dcast(seurat.object@meta.data,clonotype_id~cell.state,value.var="cell.state",length)
  # dca$Sum <- rowSums(dca[-1])
  # dca[2:(ncol(dca)-1)]<- dca[2:(ncol(dca)-1)]/dca$Sum
  # dca<-melt(dca[-ncol(dca)], id.vars = "clonotype_id")
  # dca<-dca[which(dca$clonotype_id %in% clonotype.ids),]
  desired_strings <- paste0("clonotype",clonotype.ids)
  possible_states <- unique(seurat.object$cell.state)
  temp.matrix <- as.data.frame(matrix(0,nrow=length(clonotype.ids),ncol=length(possible_states)))
  for(i in 1:nrow(temp.matrix)){
    for(j in 1:ncol(temp.matrix)){
      temp.matrix[i,j] <- length(which(seurat.object$cell.state[which(seurat.object$clonotype_id==desired_strings[i])]==possible_states[j]))/length(seurat.object$cell.state[which(seurat.object$clonotype_id==desired_strings[i])])
    }
  }
print(class(temp.matrix))
  rownames(temp.matrix) <- desired_strings
  colnames(temp.matrix) <- possible_states
  temp.matrix$clonotype_id <- desired_strings
  print(temp.matrix)

  temp.melt <- melt(temp.matrix,id.vars = "clonotype_id")
  print(temp.melt)
  stacked.ggplot<-ggplot(data=temp.melt, aes(x=clonotype_id, y=value, fill=variable)) +
    geom_bar(stat="identity")+
    ylab("Cell Counts")+
    xlab("Clonotypes")+
    labs(fill = "Cell State")
  stacked.ggplot <- stacked.ggplot + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Clonal rank") + ggplot2::theme_classic()
  return(stacked.ggplot)
}
