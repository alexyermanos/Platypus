#' Integrates VDJ and gene expression libraries by providing cluster membership seq_per_vdj object and the index of the cell in the Seurat RNA-seq object.
#' ! For platypus.version == "v3" and VDJ_GEX_matrix output the function will iterate over entries in the sample_id column of the GEX.matrix by default.
#' @param GEX.matrix For platypus.version == "v3" the GEX object from the output of the VDJ_GEX_matrix function. For platypus.version == "v2" a single seurat object from automate_GEX function after labeling cell phenotypes using the GEX_phenotype function.
#' @param clonotype.ids  For platypus.version == "v2" Output from either VDJ_analyze or VDJ_clonotype functions. This list should correspond to a single GEX.list object, in which each list element in clonotype.list is found in the GEX.object. Furthermore, these repertoires should be found in the automate_GEX library.
#' @param GEX.group.by For platypus.version == "v3". Character. Column name of the GEX.matrix@meta.data to group barplot by. Defaults to seurat_clusters
#' @param GEX.clonotypes For platypus.version == "v3". Numeric vector with ids of clonotypes to plot e.g. c(1,2,3,4). Can also be set to "topclones"
#' @param global.clonotypes Boolean. Defaults to FALSE. Set to True if clonotyping has been done across samples
#' @param platypus.version Set to either v2 or v3 automatically dependent on input format: v2 if clonotype.ids is provided and v3 if GEX.clonotypes is provided
#' @return Returns a stacked barplot that visualizes the seurat cluster membership for different cell phenotypes.
#' @export
#' @examples
#' \dontrun{
#' GEX_phenotype_per_clone_plot <- GEX_phenotype_per_clone(seurat.object = automate.gex.output[[1]], clonotype.ids= c(1,2,3,4,5))
#'}
#'
GEX_phenotype_per_clone <- function(GEX.matrix,
                                    clonotype.ids,
                                    global.clonotypes,
                                    GEX.group.by,
                                    GEX.clonotypes,
                                    platypus.version){

  value <- NULL
  variable <- NULL
  clonotype_id_10x <- NULL

  require(ggplot2)
  require(reshape2)
  require(Seurat)

  if(!missing(clonotype.ids)) platypus.version <- "v2"

  if(!missing(GEX.clonotypes)) platypus.version <- "v3"
  if(missing(global.clonotypes)) global.clonotypes <- F
  if(platypus.version == "v2"){

    seurat.object <- GEX.matrix

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

  temp.melt <- reshape2::melt(temp.matrix,id.vars = "clonotype_id")
  print(temp.melt)
  stacked.ggplot<-ggplot2::ggplot(data=temp.melt, ggplot2::aes(x=clonotype_id, y=value, fill=variable)) +
    ggplot2::geom_bar(stat="identity")+
    ggplot2::ylab("Cell Counts")+
    ggplot2::xlab("Clonotypes")+
    ggplot2::labs(fill = "Cell State")
  stacked.ggplot <- stacked.ggplot + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Clonal rank") + ggplot2::theme_classic()
  return(stacked.ggplot)

  } else if(platypus.version == "v3"){

    if(global.clonotypes){GEX.matrix@meta.data$sample_id <- 1}

    plot.out.list <- list()
    for(k in 1:length(unique(GEX.matrix@meta.data$sample_id))){

    GEX.matrix@meta.data$iscursample <- NULL
    GEX.matrix@meta.data$iscursample <- GEX.matrix@meta.data$sample_id == unique(GEX.matrix@meta.data$sample_id)[k]

    seurat.object <- subset(GEX.matrix, cells = colnames(GEX.matrix)[which(GEX.matrix$iscursample == TRUE)])

    #get clonotypes to plot
    if(GEX.clonotypes == "topclones"){ #choose 10 top expanded clones
      temp_choice <- seurat.object@meta.data[,c("clonotype_id_10x", "clonotype_frequency")]
      temp_choice <- subset(temp_choice, is.na(clonotype_id_10x) == F)
      temp_choice$clonotype_frequency <- as.numeric(temp_choice$clonotype_frequency)
      temp_choice <- temp_choice[order(temp_choice$clonotype_frequency, decreasing = T),]

      temp_choice <- temp_choice[c(1:10),]
      desired_strings <- unique(temp_choice$clonotype_id_10x)
      print(paste0("Selected ", desired_strings))
    } else { desired_strings <- paste0("clonotype",GEX.clonotypes)}
    if(any(!desired_strings %in% seurat.object$clonotype_id_10x)) stop("At least one provided clonotype does not exist in the GEX.matrix object")
    #get group ides to plot
    seurat.object$cell.state <- seurat.object@meta.data[,GEX.group.by]

    if(any(is.na(seurat.object$cell.state))){
    #subset further to exclude NA values in the grouping column
    seurat.object$groupisna <- is.na(seurat.object$cell.state) ==T
    seurat.object <- subset(seurat.object, cells = colnames(seurat.object)[which(seurat.object$groupisna == TRUE)])
    }

    possible_states <- unique(seurat.object$cell.state)
    temp.matrix <- as.data.frame(matrix(0,nrow=length(desired_strings),ncol=length(possible_states)))
    for(i in 1:nrow(temp.matrix)){
      for(j in 1:ncol(temp.matrix)){
        temp.matrix[i,j] <- length(which(seurat.object$cell.state == possible_states[j] & seurat.object$clonotype_id_10x==desired_strings[i]))/length(seurat.object$cell.state[which(seurat.object$clonotype_id_10x==desired_strings[i])]) *100
      }
    }

    rownames(temp.matrix) <- desired_strings
    colnames(temp.matrix) <- possible_states
    temp.matrix$clonotype_id <- desired_strings


    temp.melt <- reshape2::melt(temp.matrix,id.vars = "clonotype_id")

    if(is.null(levels(seurat.object$cell.state)) == F){
    temp.melt$variable <- ordered(as.factor(temp.melt$variable), levels = levels(seurat.object$cell.state))
    }

    temp.melt$clonotype_id <- ordered(as.factor(temp.melt$clonotype_id), levels = desired_strings)

    print(temp.melt)
    stacked.ggplot<- ggplot2::ggplot(data=temp.melt, ggplot2::aes(x=clonotype_id, y=value, fill=variable)) +
      ggplot2::geom_bar(stat="identity")+
      ggplot2::ylab("Cell Counts")+
      ggplot2::xlab("Clonotypes")+
      ggplot2::labs(fill = "Cell State")
    stacked.ggplot <- stacked.ggplot + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Clonal rank") + ggplot2::theme_classic() +ggplot2::scale_fill_manual(values = rainbow(length(unique(temp.melt$variable)))) + ggplot2::ggtitle(label = paste0(unique(GEX.matrix@meta.data$sample_id)[k]))
    plot.out.list[[k]] <- stacked.ggplot
    } #end loop over sample_ids
    return(plot.out.list)
  }
}



