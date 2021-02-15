#' Produces a barplot with the most frequently used IgH and IgK/L Vgenes.
#' @param clonotype.list Output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire.
#' @param HC,gene.number Numeric value indicating the top genes to be dispayed. If this number is higher than the total number of unique HC V genes in the VDJ repertoire, then this number is equal to the number of unique HC V genes.
#' @param LC.Vgene Logical indicating whether to make a barplot of the LC V genes distribution. Default is set to FALSE.
#' @param LC.gene.number Numeric value indicating the top genes to be dispayed. If this number is higher than the total number of unique LC V genes in the VDJ repertoire, then this number is equal to the number of unique LC V genes.
#' @return Returns a list of ggplot objects which show the distribution of IgH and IgK/L V genes for the most used V genes.
#' @export

VDJ_Vgene_usage_barplot <- function(clonotype.list, HC.gene.number, LC.Vgene, LC.gene.number){
  require(ggplot2)
  HC_Vgene_usage <- list()
  HC_Vgene_usage_plot <- list()
  LC_Vgene_usage <- list()
  LC_Vgene_usage_plot <- list()
  if(missing(LC.Vgene)) LC.Vgene <- FALSE
  
#  if(length(clonotype.list)>1){
    
    for (i in 1:length(clonotype.list)){
      
      HC_Vgene_usage[[i]] <- as.data.frame(table(clonotype.list[[i]]$HC_vgene))
      colnames(HC_Vgene_usage[[i]]) <- c("Vgene", "Frequency")
      for (j in 1:nrow(HC_Vgene_usage[[i]])){
        HC_Vgene_usage[[i]]$Percentage[j] <- HC_Vgene_usage[[i]]$Frequency[j]/sum(HC_Vgene_usage[[i]]$Frequency)*100
      }
      actual.gene.number <- HC.gene.number
      if(HC.gene.number > nrow(HC_Vgene_usage[[i]])){
        actual.gene.number <- nrow(HC_Vgene_usage[[i]])
      }
      #Rank based on most used V gene
      ranks <- order(-HC_Vgene_usage[[i]]$Frequency)
      HC_Vgene_usage[[i]] <- HC_Vgene_usage[[i]][ranks,]
      
      HC_Vgene_usage_plot[[i]] <- ggplot(HC_Vgene_usage[[i]][1:actual.gene.number,], aes(x=Vgene, y=Percentage, fill=Vgene)) + geom_bar(stat="identity", size=0.5, width=0.6, color="black", show.legend = FALSE) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% of unique clones") + xlab("") + theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0,0)) + scale_color_gradient(low="blue", high="red") + ggtitle(paste0("IgH V gene usage - sample ",i)) + theme(text = element_text(size = 12))
      
    }
    
#  }
  
    if(LC.Vgene==T){
      
      for (i in 1:length(clonotype.list)){
        
        LC_Vgene_usage[[i]] <- as.data.frame(table(clonotype.list[[i]]$LC_vgene))
        colnames(LC_Vgene_usage[[i]]) <- c("Vgene", "Frequency")
        for (j in 1:nrow(LC_Vgene_usage[[i]])){
          LC_Vgene_usage[[i]]$Percentage[j] <- LC_Vgene_usage[[i]]$Frequency[j]/sum(LC_Vgene_usage[[i]]$Frequency)*100
        }
        actual.gene.number <- LC.gene.number
        if(LC.gene.number > nrow(LC_Vgene_usage[[i]])){
          actual.gene.number <- nrow(LC_Vgene_usage[[i]])
        }
        #Rank based on most used V gene
        ranks <- order(-LC_Vgene_usage[[i]]$Frequency)
        LC_Vgene_usage[[i]] <- LC_Vgene_usage[[i]][ranks,]
        
        LC_Vgene_usage_plot[[i]] <- ggplot(LC_Vgene_usage[[i]][1:actual.gene.number,], aes(x=Vgene, y=Percentage, fill=Vgene)) + geom_bar(stat="identity", size=0.5, width=0.6, color="black", show.legend = FALSE) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% of unique clones") + xlab("") + theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0,0)) + scale_color_gradient(low="blue", high="red") + ggtitle(paste0("IgK/L V gene usage - sample ",i)) + theme(text = element_text(size = 12))
        
      }
      
    }
  Vgene_usage_plot <- do.call(c, list(HC_Vgene_usage_plot, LC_Vgene_usage_plot))
  
#  if(length(clonotype.list)=1)
    
    
  return(Vgene_usage_plot)
  
}


