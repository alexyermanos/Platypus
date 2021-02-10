#' Produces a stacked barplot with the fraction of the most frequently used IgH and IgK/L Vgenes. This function can be used in combination with the VDJ_Vgene_usage_barplot to vizualize V gene usage per sample and among samples.
#' @param clonotype.list Output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire.
#' @param HC.gene.number Numeric value indicating the top genes to be dispayed. If this number is higher than the total number of unique HC V genes in the VDJ repertoire, then this number is equal to the number of unique HC V genes.
#' @param Fraction.HC Numeric value indicating the minimum fraction of clones expressing a particular HC V gene. If the usage of a particular gene is below this value, then this gene is excluded. If the usage of a particular gene is above this value even in one sample, then this gene is included in the analysis. Default value is set to 0, thus all genes are selected.
#' @param LC.Vgene Logical indicating whether to make a barplot of the LC V gene distribution. Default is set to FALSE.
#' @param LC.gene.number Numeric value indicating the top genes to be dispayed. If this number is higher than the total number of unique LC V genes in the VDJ repertoire, then this number is equal to the number of unique LC V genes.
#' @param Fraction.LC Numeric value indicating the minimum fraction of clones expressing a particular LC V gene. If the usage of a particular gene is below this value, then this gene is excluded. If the usage of a particular gene is above this value even in one sample, then this gene is included in the analysis. Default value is set to 0, thus all genes are selected.
#' @return Returns a list of ggplot objects which show the stacked distribution of IgH and IgK/L V genes for the most used V genes. Returns an empty plot if the Fraction.HC or Fraction.LC that were selected were too high, resulting in the exclusion of all the genes.
#' \dontrun{
#' example.vdj.vgene_usage <- VDJ_Vgene_usage_barplot_stacked(clonotype.list = covid_vdj_repertoire_bcells, LC.Vgene = T,HC.gene.number = 15, Fraction.HC = 1)
#'}
#' @export
VDJ_Vgene_usage_stacked_barplot <- function(clonotype.list, HC.gene.number, Fraction.HC, LC.Vgene, LC.gene.number, Fraction.LC){
  require(ggplot2)
  if (missing(Fraction.HC)) Fraction.HC <- 0
  if (missing(HC.gene.number)) HC.gene.number <- 10
  if (missing(LC.Vgene)) LC.Vgene <- FALSE
  if (missing(Fraction.LC)) Fraction.LC <- 0
  if (missing(LC.gene.number)) LC.gene.number <- 10

  Vgene_usage_plot <- list()
  HC_Vgene_usage <- list()
  HC_Vgene_usage_top <- list()
  HC_Vgene_usage_fraction <- list()

  for (i in 1:length(clonotype.list)){
    #Calculate the V gene usage frequencies for each sample on the clonal lebel
    HC_Vgene_usage[[i]] <- as.data.frame(table(clonotype.list[[i]]$HC_vgene))
    colnames(HC_Vgene_usage[[i]]) <- c("Vgene", "Frequency")
    for (j in 1:nrow(HC_Vgene_usage[[i]])){
      HC_Vgene_usage[[i]]$Percentage[j] <- HC_Vgene_usage[[i]]$Frequency[j]/sum(HC_Vgene_usage[[i]]$Frequency)*100
    }

    #Rank based on most used V gene per sample
    ranks <- order(-HC_Vgene_usage[[i]]$Frequency)
    HC_Vgene_usage[[i]] <- HC_Vgene_usage[[i]][ranks,]
    #Find the genes which are used more than Fraction.HC of the time
    HC_Vgene_usage_fraction[[i]] <- HC_Vgene_usage[[i]][which(HC_Vgene_usage[[i]]$Percentage > Fraction.HC),]
  }
  #Find the union of unique gene names between all the samples
  HC_Vgene_names <- do.call("rbind",HC_Vgene_usage_fraction)
  HC_Vgene_names <- as.vector(unique(HC_Vgene_names$Vgene))

  #Filter out all those genes that are used less than the Fraction.HC in all samples
  for(i in 1:length(clonotype.list)){
    HC_Vgene_usage[[i]] <- HC_Vgene_usage[[i]][HC_Vgene_usage[[i]]$Vgene %in% HC_Vgene_names,]
  }

  #Determine the mean gene usage between the samples

  HC_Vgene_usage_all <- do.call("rbind", HC_Vgene_usage)
  unique_Vgene <- unique(HC_Vgene_usage_all$Vgene)

  HC_Vgene_usage_all_unique <- data.frame("Vgene"= rep("", length(unique_Vgene)), "Frequency"= rep(0, length(unique_Vgene)))
  HC_Vgene_usage_all_unique$Vgene <- unique_Vgene

  for (i in 1:length(unique_Vgene)){
    HC_Vgene_usage_all_unique[i,"Percentage"] <- mean(HC_Vgene_usage_all[which(HC_Vgene_usage_all$Vgene == HC_Vgene_usage_all_unique$Vgene[i]),]$Percentage)
  }
  # Rank the V genes based on the most used ones
  ranks <- order(-HC_Vgene_usage_all_unique$Frequency)
  HC_Vgene_usage_all_unique <- HC_Vgene_usage_all_unique[ranks,]
  top_Vgenes <- HC_Vgene_usage_all_unique$Vgene[1:HC.gene.number]

  for(i in 1:length(clonotype.list)){
    HC_Vgene_usage_top[[i]] <- data.frame("Vgene"=rep("", length(top_Vgenes)), "Frequency"=rep(0, length(top_Vgenes)))
    HC_Vgene_usage_top[[i]]$Vgene <- top_Vgenes
    for (j in 1:nrow(HC_Vgene_usage_top[[i]])){
      if (HC_Vgene_usage_top[[i]]$Vgene[j] %in% HC_Vgene_usage[[i]]$Vgene == TRUE){
        HC_Vgene_usage_top[[i]]$Frequency[j] <- HC_Vgene_usage[[i]][which(HC_Vgene_usage[[i]]$Vgene == paste0(HC_Vgene_usage_top[[i]]$Vgene[j])),]$Frequency
      }
    }
    for (j in 1:nrow(HC_Vgene_usage_top[[i]])){
      HC_Vgene_usage_top[[i]]$Percentage[j] <- HC_Vgene_usage_top[[i]]$Frequency[j]/sum(HC_Vgene_usage_top[[i]]$Frequency)*100
    }
  }

  #Assign the sample id
  for (i in 1:length(HC_Vgene_usage_top)){
    HC_Vgene_usage_top[[i]]$Sample <- as.character(i)
  }

  #And then bind together for plotting
  plotting_df <- do.call("rbind",HC_Vgene_usage_top)
  #Take out the NAs introduced because of the HC.gene.number and the Fraction.HC.
  plotting_df <- plotting_df[!is.na(plotting_df$Vgene), ]

  Vgene_usage_plot[[1]] <- ggplot(plotting_df, aes(fill = Vgene, y=Frequency, x=Sample)) +
    geom_bar(position="fill", stat="identity", color="black", width = 0.7) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("% of unique clones") + scale_y_continuous(expand = c(0,0))

  if(LC.Vgene ==TRUE){

    LC_Vgene_usage <- list()
    LC_Vgene_usage_top <- list()
    LC_Vgene_usage_fraction <- list()

    for (i in 1:length(clonotype.list)){

      LC_Vgene_usage[[i]] <- as.data.frame(table(clonotype.list[[i]]$LC_vgene))
      colnames(LC_Vgene_usage[[i]]) <- c("Vgene", "Frequency")
      for (j in 1:nrow(LC_Vgene_usage[[i]])){
        LC_Vgene_usage[[i]]$Percentage[j] <- LC_Vgene_usage[[i]]$Frequency[j]/sum(LC_Vgene_usage[[i]]$Frequency)*100
      }

      #Rank based on most used V gene
      ranks <- order(-LC_Vgene_usage[[i]]$Frequency)
      LC_Vgene_usage[[i]] <- LC_Vgene_usage[[i]][ranks,]

      LC_Vgene_usage_fraction[[i]] <- LC_Vgene_usage[[i]][which(LC_Vgene_usage[[i]]$Percentage > Fraction.LC),]

    }

    #Find the union of unique gene names between all the samples
    LC_Vgene_names <- do.call("rbind",LC_Vgene_usage_fraction)
    LC_Vgene_names <- as.vector(unique(LC_Vgene_names$Vgene))

    #Filter out all those genes that are used less than the Fraction.HC in all samples
    for(i in 1:length(clonotype.list)){
      LC_Vgene_usage[[i]] <- LC_Vgene_usage[[i]][LC_Vgene_usage[[i]]$Vgene %in% LC_Vgene_names,]
    }

    LC_Vgene_usage_all <- do.call("rbind", LC_Vgene_usage)
    unique_Vgene <- unique(LC_Vgene_usage_all$Vgene)

    LC_Vgene_usage_all_unique <- data.frame("Vgene"= rep("", length(unique_Vgene)), "Frequency"= rep(0, length(unique_Vgene)))
    LC_Vgene_usage_all_unique$Vgene <- unique_Vgene

    for (i in 1:length(unique_Vgene)){
      LC_Vgene_usage_all_unique[i,"Percentage"] <- mean(LC_Vgene_usage_all[which(LC_Vgene_usage_all$Vgene == LC_Vgene_usage_all_unique$Vgene[i]),]$Percentage)
    }
    # Rank the V genes based on the most used ones
    ranks <- order(-LC_Vgene_usage_all_unique$Frequency)
    LC_Vgene_usage_all_unique <- LC_Vgene_usage_all_unique[ranks,]
    top_Vgenes <- LC_Vgene_usage_all_unique$Vgene[1:LC.gene.number]

    for(i in 1:length(clonotype.list)){
      LC_Vgene_usage_top[[i]] <- data.frame("Vgene"=rep("", length(top_Vgenes)), "Frequency"=rep(0, length(top_Vgenes)))
      LC_Vgene_usage_top[[i]]$Vgene <- top_Vgenes
      for (j in 1:nrow(LC_Vgene_usage_top[[i]])){
        if (LC_Vgene_usage_top[[i]]$Vgene[j] %in% LC_Vgene_usage[[i]]$Vgene == TRUE){
          LC_Vgene_usage_top[[i]]$Frequency[j] <- LC_Vgene_usage[[i]][which(LC_Vgene_usage[[i]]$Vgene == paste0(LC_Vgene_usage_top[[i]]$Vgene[j])),]$Frequency
        }
      }
      for (j in 1:nrow(LC_Vgene_usage_top[[i]])){
        LC_Vgene_usage_top[[i]]$Percentage[j] <- LC_Vgene_usage_top[[i]]$Frequency[j]/sum(LC_Vgene_usage_top[[i]]$Frequency)*100
      }
    }

    #Assign the sample id
    for (i in 1:length(LC_Vgene_usage_top)){
      LC_Vgene_usage_top[[i]]$Sample <- as.character(i)
    }

    #And then bind together
    plotting_df <- do.call("rbind",LC_Vgene_usage_top)
    plotting_df <- plotting_df[!is.na(plotting_df$Vgene), ]

    Vgene_usage_plot[[2]] <-ggplot(plotting_df, aes(fill = Vgene, y=Frequency, x=Sample)) +
      geom_bar(position="fill", stat="identity", color="black", width = 0.7) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ylab("% of unique clones") + scale_y_continuous(expand = c(0,0))

  }

  return(Vgene_usage_plot)
}