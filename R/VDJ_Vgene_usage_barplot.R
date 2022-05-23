#'V(D)J gene usage barplots
#'
#' @description Produces a barplot with the most frequently used IgH and IgK/L Vgenes.
#' @param VDJ Either (for platypus version "v2") output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire, OR (for platypus version "v3") the the VDJ matrix output of the VDJ_GEX_matrix() function (VDJ.GEX.matrix.output[[1]])
#' @param group.by Character. Defaults to "sample_id". Column name of VDJ to group plot by.
#' @param HC.gene.number Numeric value indicating the top genes to be dispayed. If this number is higher than the total number of unique HC V genes in the VDJ repertoire, then this number is equal to the number of unique HC V genes.
#' @param LC.Vgene Logical indicating whether to make a barplot of the LC V genes distribution. Default is set to FALSE.
#' @param LC.gene.number Numeric value indicating the top genes to be dispayed. If this number is higher than the total number of unique LC V genes in the VDJ repertoire, then this number is equal to the number of unique LC V genes.
#' @param platypus.version Character. Defaults to "v3". Can be "v2" or "v3" dependent on the input format
#' @return Returns a list of ggplot objects which show the distribution of IgH and IgK/L V genes for the most used V genes.
#' @param is.bulk logical value indicating whether the VDJ input was generated from bulk-sequencing data using the bulk_to_vgm function. If is.bulk = T, the VDJ_Vgene_usage_barplot function is compatible for use with bulk data. Defaults to False (F).
#' @export
#' @examples
#' \dontrun{
#' VDJ_Vgene_usage_barplot(VDJ = Platypus::small_vgm[[1]],
#' HC.gene.number = 2, platypus.version = "v3")
#'}

VDJ_Vgene_usage_barplot <- function(VDJ,
                                    group.by,
                                    HC.gene.number,
                                    LC.Vgene,
                                    LC.gene.number,
                                    platypus.version,
                                    is.bulk){

  Vgene <- NULL
  Percentage <- NULL
  Nr_of_VDJ_chains <- NULL
  Nr_of_VJ_chains <- NULL
  sample_id <- NULL
  HC_Vgene_usage <- list()
  HC_Vgene_usage_plot <- list()
  LC_Vgene_usage <- list()
  LC_Vgene_usage_plot <- list()

  if(missing(is.bulk)) is.bulk <- F
  if(missing(LC.Vgene)) LC.Vgene <- FALSE
  if(missing(platypus.version)) platypus.version <- "v3"
  if(missing(HC.gene.number)) HC.gene.number <- 10
  #naming compatibility
  VDJ.matrix <- VDJ
  VDJ <- NULL
  if(platypus.version == "v2"){

    clonotype.list <- VDJ.matrix

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

      HC_Vgene_usage_plot[[i]] <- ggplot2::ggplot(HC_Vgene_usage[[i]][1:actual.gene.number,], ggplot2::aes(x=Vgene, y=Percentage, fill=Vgene)) + ggplot2::geom_bar(stat="identity", size=0.5, width=0.6, color="black", show.legend = FALSE) + cowplot::theme_cowplot() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + ggplot2::ylab("% of unique clones") + ggplot2::xlab("") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::scale_color_gradient(low="blue", high="red") + ggplot2::ggtitle(paste0("sample ",i)) + ggplot2::theme(text = ggplot2::element_text(size = 12))

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

        LC_Vgene_usage_plot[[i]] <- ggplot2::ggplot(LC_Vgene_usage[[i]][1:actual.gene.number,], ggplot2::aes(x=Vgene, y=Percentage, fill=Vgene)) + ggplot2::geom_bar(stat="identity", size=0.5, width=0.6, color="black", show.legend = FALSE) + cowplot::theme_cowplot() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + ggplot2::ylab("% of unique clones") + ggplot2::xlab("") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::scale_color_gradient(low="blue", high="red") + ggplot2::ggtitle(paste0("sample ",i)) + ggplot2::theme(text = ggplot2::element_text(size = 12))

      }

    }
  Vgene_usage_plot <- do.call(c, list(HC_Vgene_usage_plot, LC_Vgene_usage_plot))



  return(Vgene_usage_plot)

  } else if(platypus.version == "v3"){

    if(is.bulk == F){
      #filtering for max 1VDJ 1VJ chain
      VDJ.matrix <- subset(VDJ.matrix, Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)
    }

    if(missing(group.by)) group.by <- "sample_id"
    if(group.by != "sample_id"){
      if(group.by %in% names(VDJ.matrix)){
        VDJ.matrix$sample_id <- as.character(VDJ.matrix[,group.by])
        if(any(is.na(VDJ.matrix$sample_id)) == T){
          VDJ.matrix <- VDJ.matrix[!is.na(VDJ.matrix$sample_id),]
          warning(paste0("Filtered out cells with 'NA' in grouping column"))
        }
        message(paste0("Grouping by: ", group.by))
      } else {
        warning(paste0("Group_id '",group.by, "' was not found in VDJ. Grouping by 'sample_id'"))}
    }

    clonotype.list <- list()
    for(i in 1:length(unique(VDJ.matrix$sample_id))){
      clonotype.list[[i]] <- subset(VDJ.matrix, sample_id == unique(VDJ.matrix$sample_id)[i])
      #removing extra cells cells to leave only 1 per clonotype
      clonotype.list[[i]] <- clonotype.list[[i]][duplicated(clonotype.list[[i]]$clonotype_id_10x) == F,]
    }
    names(clonotype.list) <- unique(VDJ.matrix$sample_id)
    message(paste0("Sample order: ", paste0(unique(VDJ.matrix$sample_id), collapse = " ; ")))


    for (i in 1:length(clonotype.list)){

      HC_Vgene_usage[[i]] <- as.data.frame(table(clonotype.list[[i]]$VDJ_vgene))
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

      HC_Vgene_usage_plot[[i]] <- ggplot2::ggplot(HC_Vgene_usage[[i]][1:actual.gene.number,], ggplot2::aes(x=Vgene, y=Percentage, fill=Vgene)) + ggplot2::geom_bar(stat="identity", size=0.5, width=0.6, color="black", show.legend = FALSE) + cowplot::theme_cowplot() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + ggplot2::ylab("% of unique clones") + ggplot2::xlab("") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::scale_color_gradient(low="blue", high="red") + ggplot2::ggtitle(paste0("sample ",names(clonotype.list)[i])) + ggplot2::theme(text = ggplot2::element_text(size = 12))

    }

    #  }

    if(LC.Vgene==T){

      for (i in 1:length(clonotype.list)){

        LC_Vgene_usage[[i]] <- as.data.frame(table(clonotype.list[[i]]$VJ_vgene))
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

        LC_Vgene_usage_plot[[i]] <- ggplot2::ggplot(LC_Vgene_usage[[i]][1:actual.gene.number,], ggplot2::aes(x=Vgene, y=Percentage, fill=Vgene)) + ggplot2::geom_bar(stat="identity", size=0.5, width=0.6, color="black", show.legend = FALSE) + cowplot::theme_cowplot() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + ggplot2::ylab("% of unique clones") + ggplot2::xlab("") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::scale_color_gradient(low="blue", high="red") + ggplot2::ggtitle(paste0("sample ", names(clonotype.list)[i])) + ggplot2::theme(text = ggplot2::element_text(size = 12))

      }

    }
    Vgene_usage_plot <- do.call(c, list(HC_Vgene_usage_plot, LC_Vgene_usage_plot))

    return(Vgene_usage_plot)
  }
}

