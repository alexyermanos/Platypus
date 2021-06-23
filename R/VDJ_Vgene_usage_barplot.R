#' Produces a barplot with the most frequently used IgH and IgK/L Vgenes.
#' @param VDJ.matrix Either (for platypus version "v2") output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire, OR (for platypus version "v3") the the VDJ matrix output of the VDJ_GEX_matrix() function (normally VDJ.GEX.matrix.output[[1]])
#' @param HC.gene.number Numeric value indicating the top genes to be dispayed. If this number is higher than the total number of unique HC V genes in the VDJ repertoire, then this number is equal to the number of unique HC V genes.
#' @param LC.Vgene Logical indicating whether to make a barplot of the LC V genes distribution. Default is set to FALSE.
#' @param LC.gene.number Numeric value indicating the top genes to be dispayed. If this number is higher than the total number of unique LC V genes in the VDJ repertoire, then this number is equal to the number of unique LC V genes.
#' @param platypus.version Character. Defaults to "v2". Can be "v2" or "v3" dependent on the input format
#' @return Returns a list of ggplot objects which show the distribution of IgH and IgK/L V genes for the most used V genes.
#' @export

VDJ_Vgene_usage_barplot <- function(VDJ.matrix, HC.gene.number, LC.Vgene, LC.gene.number, platypus.version){

  Vgene <- NULL
  Percentage <- NULL
  Nr_of_VDJ_chains <- NULL
  Nr_of_VJ_chains <- NULL
  sample_id <- NULL

  require(ggplot2)
  HC_Vgene_usage <- list()
  HC_Vgene_usage_plot <- list()
  LC_Vgene_usage <- list()
  LC_Vgene_usage_plot <- list()
  if(missing(LC.Vgene)) LC.Vgene <- FALSE
  if(missing(platypus.version)) platypus.version <- "v2"

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

      HC_Vgene_usage_plot[[i]] <- ggplot2::ggplot(HC_Vgene_usage[[i]][1:actual.gene.number,], ggplot2::aes(x=Vgene, y=Percentage, fill=Vgene)) + ggplot2::geom_bar(stat="identity", size=0.5, width=0.6, color="black", show.legend = FALSE) + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + ggplot2::ylab("% of unique clones") + ggplot2::xlab("") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::scale_color_gradient(low="blue", high="red") + ggplot2::ggtitle(paste0("IgH V gene usage - sample ",i)) + ggplot2::theme(text = ggplot2::element_text(size = 12))

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

        LC_Vgene_usage_plot[[i]] <- ggplot2::ggplot(LC_Vgene_usage[[i]][1:actual.gene.number,], ggplot2::aes(x=Vgene, y=Percentage, fill=Vgene)) + ggplot2::geom_bar(stat="identity", size=0.5, width=0.6, color="black", show.legend = FALSE) + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + ggplot2::ylab("% of unique clones") + ggplot2::xlab("") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::scale_color_gradient(low="blue", high="red") + ggplot2::ggtitle(paste0("IgK/L V gene usage - sample ",i)) + ggplot2::theme(text = ggplot2::element_text(size = 12))

      }

    }
  Vgene_usage_plot <- do.call(c, list(HC_Vgene_usage_plot, LC_Vgene_usage_plot))



  return(Vgene_usage_plot)

  } else if(platypus.version == "v3"){

    #filtering for max 1VDJ 1VJ chain
    VDJ.matrix <- subset(VDJ.matrix, Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)


    clonotype.list <- list()
    for(i in 1:length(unique(VDJ.matrix$sample_id))){
      clonotype.list[[i]] <- subset(VDJ.matrix, sample_id == unique(VDJ.matrix$sample_id)[i])
      #removing extra cells cells to leave only 1 per clonotype
      clonotype.list[[i]] <- clonotype.list[[i]][duplicated(clonotype.list[[i]]$clonotype_id_10x) == F,]
    }
    names(clonotype.list) <- unique(VDJ.matrix$sample_id)
    print(paste0("Sample order: ", paste0(unique(VDJ.matrix$sample_id), collapse = " ; ")))


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

      HC_Vgene_usage_plot[[i]] <- ggplot2::ggplot(HC_Vgene_usage[[i]][1:actual.gene.number,], ggplot2::aes(x=Vgene, y=Percentage, fill=Vgene)) + ggplot2::geom_bar(stat="identity", size=0.5, width=0.6, color="black", show.legend = FALSE) + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + ggplot2::ylab("% of unique clones") + ggplot2::xlab("") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::scale_color_gradient(low="blue", high="red") + ggplot2::ggtitle(paste0("IgH V gene usage - sample ",names(clonotype.list)[i])) + ggplot2::theme(text = ggplot2::element_text(size = 12))

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

        LC_Vgene_usage_plot[[i]] <- ggplot2::ggplot(LC_Vgene_usage[[i]][1:actual.gene.number,], ggplot2::aes(x=Vgene, y=Percentage, fill=Vgene)) + ggplot2::geom_bar(stat="identity", size=0.5, width=0.6, color="black", show.legend = FALSE) + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + ggplot2::ylab("% of unique clones") + ggplot2::xlab("") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::scale_color_gradient(low="blue", high="red") + ggplot2::ggtitle(paste0("IgK/L V gene usage - sample ", names(clonotype.list)[i])) + ggplot2::theme(text = ggplot2::element_text(size = 12))

      }

    }
    Vgene_usage_plot <- do.call(c, list(HC_Vgene_usage_plot, LC_Vgene_usage_plot))

    return(Vgene_usage_plot)


  }
}
