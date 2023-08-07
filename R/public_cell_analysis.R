#'public_cell_analysis function:
#'Creates Plots using the find_public_cells() output vgm. Which type of plot and analysis should be done can be chosen by selecting the desired analysis.type input. Default is set to "overlap type counts"
#'
#'1. Plot: Overlap Type Counts per Sample (CDRH3, CDRL3 and CDRL3.CDRH3 overlaps) - analysis.type= "overlap type counts"
#'
#'2. Plot: Public Cell Numbers - analysis.type = "public cell numbers"
#'Note: It is likely that a cell which has a CDRH3 overlap has a CDRL3 overlap with a cell in the public data as well, but is not listed as a CDRL3_CDRH3_aa overlap. Because cells are only assigned a complete overlap if both (heavy and light chain) overlaps refer to the same cell in the public data.
#'Some CDRL3 overlaps might therefore be overwritten with CDRH3 overlaps due to their higher significance.
#'
#'3. Plot: Public Clone Numbers - analysis.type = "public clone numbers"
#'
#'4. Plot: Clonal expansion colored by overlap type - analysis.type = "clonal expansion by overlap type"
#'Note: Here you can use an additional input called: include.light.chain
#'If include.light.chain = TRUE, also light chain overlaps are colored in the clonal expansion plots. Default is set to FALSE
#'Also, this function returns a pdf of clonal expansion plots (all samples) instead of only one or two plots, if save.plot = TRUE
#'
#'5. Isotype usage of public cells - analysis.type = "isotype usage of public cells"
#'Creates stacked barplots of the isotype usage of all public cells
#'
#'6. Transcriptional cluster usage of public cells - analysis.type = "transcriptional cluster usage of public cells"
#'Creates stacked barplots of the transcriptional cluster usage of public cells
#'
#'@param find_public_cells_output Object, output of find_public_cells function
#'@param analysis.type Character, defines what analysis should be performed. Options are: "overlap type counts", "public cell numbers", "public clone numbers",  "clonal expansion by overlap type", "isotype usage of public cells" and "transcriptional cluster usage of public cells". Default is "overlap type counts".
#'@param inlcude.light.chain Logical, default is FALSE. Only needed if analysis.type = "clonal expansion by overlap type". If set to TRUE, also cdr3_aa light chain overlaps are colored in the created clonal expansion plots.
#'@param save.pot Logical, default is FALSE. If set to TRUE, created plots are saved in the current working directory.
#'



public_cell_analysis <- function(normalize_y_n,find_public_cells_output,
                                 analysis.type,
                                 include.light.chain,
                                 save.plot,
                                 nr_clones){

  if(missing(analysis.type)) analysis.type <- "overlap type counts"
  if(missing(include.light.chain)) include.light.chain <- FALSE
  if(missing(save.plot)) save.plot <- FALSE
  if(missing(nr_clones)) nr_clones <- '30'



  #To create plots load packages
  #library(viridis)
  #library(hrbrthemes)

  vgm <- find_public_cells_output[[1]]
  vdj <- vgm[[1]]

  public_clones <- find_public_cells_output[[2]]


  sample_ids <- unique(vdj$sample_id)

  # 1. Overlap Type Counts per Sample (CDRH3, CDRL3 and CDRL3.CDRH3 overlaps)
  if(analysis.type == "overlap type counts"){

    overlap_type_counts <- NULL

    for(i in 1:length(sample_ids)){
      df <- NULL
      sample <- subset(vdj, vdj$sample_id == sample_ids[[i]])
      df <- as.data.frame(table(sample$overlap_type))
      df$scaled_counts <- round((df$Freq)/(nrow(sample))*100, digits = 2)
      df$sample_id <- sample_ids[[i]]

      names(df) <- c("Overlap Type", "counts", "scaled_counts", "sample")

      overlap_type_counts <- rbind(overlap_type_counts, df)


    }
    #Turn Overlap Type column into factor for better visualization
    overlap_type_counts$`Overlap Type` <- as.factor(overlap_type_counts$`Overlap Type`)
    overlap_type_counts$`Overlap Type` <- factor(overlap_type_counts$`Overlap Type`, levels = c("No overlap", "CDRL3_aa", "CDRH3_aa", "CDRL3_CDRH3_aa"))


    #Create a plot of the public cell counts split by sample id
    overlap_types_barplot <- ggplot(overlap_type_counts, ggplot2::aes(fill = `Overlap Type`, y = counts, x = `sample`)) + ggplot2::geom_bar(position="stack", stat="identity") +ggtitle("Overlap Type Counts - Single Cell Level")+ xlab("") + ylab("") + viridis::scale_fill_viridis(discrete = T) + theme_ipsum() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))  + ggplot2::geom_text(aes(label=counts), position = position_stack(vjust= 0.8),colour = "red", size = 4) #Add cell numbers

    #Create a plot of the public cell frequency split by sample id
    overlap_types_scaled_barplot <- ggplot(overlap_type_counts, ggplot2::aes(fill = `Overlap Type`, y = scaled_counts, x = `sample`)) + ggplot2::geom_bar(position="stack", stat="identity") +ggtitle("Overlap Type Frequency - Single Cell Level")+ xlab("") + ylab("") + viridis::scale_fill_viridis(discrete = T) + theme_ipsum() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::geom_text(aes(label=paste0(scaled_counts, "%")), position = ggplot2::position_stack(vjust= 0.8),colour = "red", size = 4) #

    if(save.plot == TRUE){
      ggplot2::ggsave(overlap_types_barplot, filename = "overlap_types_barplot.png", dpi = 400, width = 12.0, height = 9.48)
      ggplot2::ggsave(overlap_types_scaled_barplot, filename = "overlap_types_scaled_barplot.png", dpi = 400, width = 12.0, height = 9.48)
    }

    tmp <- list(overlap_types_barplot, overlap_types_scaled_barplot)
    names(tmp) <- c("overlap_types_barplot", "overlap_types_scaled_barplot")

    return(tmp)
  }


  # 2. Public Cell Numbers
  if(analysis.type == "public cell numbers"){

    #Create a dataframe containing sample id information and public cell counts (scaled and unscaled) for plotting
    public_cell_numbers <- NULL

    for(i in 1:length(sample_ids)){
      df <- NULL
      sample <- subset(vdj, vdj$sample_id == sample_ids[[i]])
      df <- as.data.frame(table(sample$public_clonotype))
      df$scaled_counts <- round((df$Freq)/(nrow(sample))*100, digits = 2)
      df$sample_id <- sample_ids[[i]]

      names(df) <- c("Public", "counts", "scaled_counts", "sample")

      public_cell_numbers <- rbind(public_cell_numbers, df)
    }

    #Create a plot of the public cell counts split by sample id
    public_cells_barplot <- ggplot(public_cell_numbers, aes(fill = `Public`, y = counts, x = `sample`)) + geom_bar(position="stack", stat="identity") +ggtitle("Public Cell Count")+ xlab("") + ylab("") + scale_fill_viridis(discrete = T) + theme_ipsum() + theme(plot.title = element_text(hjust = 0.5))  + geom_text(aes(label=counts), position = position_stack(vjust= 0.8),colour = "red", size = 5) #Add cell numbers

    #Create a plot of the public cell frequency split by sample id
    public_cells_scaled_barplot <- ggplot(public_cell_numbers, aes(fill = `Public`, y = scaled_counts, x = `sample`)) + geom_bar(position="stack", stat="identity") +ggtitle("Public Cell Frequency")+ xlab("") + ylab("") + scale_fill_viridis(discrete = T) + theme_ipsum() + theme(plot.title = element_text(hjust = 0.5)) + geom_text(aes(label=paste0(scaled_counts, "%")), position = position_stack(vjust= 0.8),colour = "red", size = 5) #


    if(save.plot == TRUE){
      ggsave(public_cells_barplot, filename = "public_cells_barplot.png", dpi = 400, width = 12.0, height = 9.48)
      ggsave(public_cells_scaled_barplot, filename = "public_cells_scaled_barplot.png", dpi = 400, width = 12.0, height = 9.48)
    }

    tmp <- list(public_cells_barplot, public_cells_scaled_barplot)
    names(tmp) <- c("public_cells_barplot", "public_cells_scaled_barplot")

    return(tmp)


    #Note: It is likely that a cell which has a CDRH3 overlap has a CDRL3 overlap with a cell in the public data as well, but is not listed as a CDRL3_CDRH3_aa overlap. Because cells are only assigned a complete overlap if both (heavy and light chain) overlaps refer to the same cell in the public data.
    #Some CDRL3 overlaps might therefore be overwritten with CDRH3 overlaps due to their higher significance.

  }


  # 3. Public Clone Numbers
  if(analysis.type == "public clone numbers"){

    public_clone_numbers <- NULL

    for(i in 1:length(sample_ids)){
      df <- NULL
      sample <- subset(vdj, vdj$sample_id == sample_ids[[i]])

      df <- as.data.frame(table(sample$public_clonotype))
      df$scaled_counts <- round((df$Freq)/(nrow(sample))*100, digits = 2)
      df$sample_id <- sample_ids[[i]]

      names(df) <- c("Public", "counts", "scaled_counts", "sample")

      public_clone_numbers <- rbind(public_clone_numbers, df)
    }

    #Create a plot of the public clone counts split by sample id
    public_clones_barplot <- ggplot(public_clone_numbers, aes(fill = `Public`, y = counts, x = `sample`)) + geom_bar(position="stack", stat="identity") +ggtitle("Public Clone Count")+ xlab("") + ylab("") + scale_fill_viridis(discrete = T) + theme_ipsum() + theme(plot.title = element_text(hjust = 0.5))  + geom_text(aes(label=counts), position = position_stack(vjust= 0.8),colour = "red", size = 5) #Add cell numbers

    #Create a plot of the public clone frequency split by sample id
    public_clones_scaled_barplot <- ggplot(public_clone_numbers, aes(fill = `Public`, y = scaled_counts, x = `sample`)) + geom_bar(position="stack", stat="identity") +ggtitle("Public Clone Frequency")+ xlab("") + ylab("") + scale_fill_viridis(discrete = T) + theme_ipsum() + theme(plot.title = element_text(hjust = 0.5)) + geom_text(aes(label=paste0(scaled_counts, "%")), position = position_stack(vjust= 0.8),colour = "red", size = 5) #


    if(save.plot == TRUE){
      ggsave(public_clones_barplot, filename = "public_clones_barplot.png", dpi = 400, width = 12.0, height = 9.48)
      ggsave(public_clones_scaled_barplot, filename = "public_clones_scaled_barplot.png", dpi = 400, width = 12.0, height = 9.48)
    }

    tmp <- list(public_clones_barplot, public_clones_scaled_barplot)
    names(tmp) <- c("public_clones_barplot", "public_clones_scaled_barplot")

    return(tmp)

  }



  # 4. Plot: Clonal expansion colored by overlap type
  if(analysis.type == "clonal expansion by overlap type"){

    vdj$overlap_type <- as.factor(vdj$overlap_type)
    vdj$overlap_type <- factor(vdj$overlap_type, levels = c("No overlap", "CDRL3_aa", "CDRH3_aa", "CDRL3_CDRH3_aa"))

    if (normalize_y_n == 'y') {


      if(include.light.chain == TRUE){
        clonal_expansion_by_overlap_type <-  VDJ_clonal_expansion_normalized_MV(VDJ = vdj,  group.by = "sample_id", color.by = "overlap_type", clones = nr_clones, isotypes.to.plot = "all", species = "Mouse", treat.incomplete.clones = "exclude", treat.incomplete.cells = "proportional", normalization_method = "split_capture", groups_of_capture = list(capture1, capture2), until_which_clone_normalize = clones_with_two_cells)
      }else{
        vdj$overlap_type[vdj$overlap_type == "CDRL3_aa"] <- "No overlap" #overwrite light chain overlap with no overlap (e.g. usefull if you have a lot of light chain overlaps)
        clonal_expansion_by_overlap_type <-  VDJ_clonal_expansion_normalized_MV(VDJ = vdj,  group.by = "sample_id", color.by = "overlap_type", clones = nr_clones, isotypes.to.plot = "all", species = "Mouse", treat.incomplete.clones = "exclude", treat.incomplete.cells = "proportional", normalization_method = "split_capture", groups_of_capture = list(capture1, capture2), until_which_clone_normalize = clones_with_two_cells)
      }

    }

    else if (normalize_y_n == 'n') {


      if(include.light.chain == TRUE){
        clonal_expansion_by_overlap_type <-  VDJ_clonal_expansion(VDJ = vdj, color.by = "overlap_type", clones = nr_clones)
      }else{
        vdj$overlap_type[vdj$overlap_type == "CDRL3_aa"] <- "No overlap" #overwrite light chain overlap with no overlap (e.g. usefull if you have a lot of light chain overlaps)
        clonal_expansion_by_overlap_type <-  VDJ_clonal_expansion(VDJ = vdj, color.by = "overlap_type", clones = nr_clones)
      }

    }

    #Change plot titles
    for(i in 1:length(clonal_expansion_by_overlap_type[[1]])){
      clonal_expansion_by_overlap_type[[1]][[i]][["labels"]][["title"]] <- paste(clonal_expansion_by_overlap_type[[1]][[i]][["labels"]][["title"]], "Clonal expansion colored by cdr3 (aa) overlap type", sep = ": ")
    }

    if(save.plot == TRUE){
      ggsave(
        filename = "clonal_expansion_by_cdr3_aa_overlap.pdf",
        plot = marrangeGrob(clonal_expansion_by_overlap_type[[1]], nrow=1, ncol=1),
        width = 15, height = 9
      )
    }

    return(clonal_expansion_by_overlap_type)
  }


  # 5. Isotype usage of public cells
  if(analysis.type == "isotype usage of public cells"){

    #Create subset containing only cells with a public clonotype
    public_cells <- subset(vdj, vdj$public == TRUE)

    #create stacked barplots per isotype usage

    isotype_usage <- NULL

    for(i in 1:length(sample_ids)){
      curr_sample <- subset(public_cells, public_cells$sample_id == sample_ids[i])

      df <- as.data.frame(table(curr_sample$VDJ_cgene))
      df$scaled_freq <- df$Freq/nrow(curr_sample)*100
      df$sample <- sample_ids[i]

      names(df) <- c("Isotype", "Counts", "Scaled_Counts", "Sample")

      isotype_usage <- rbind(isotype_usage,df)

    }



    isotype_usage_public_cells <- ggplot(isotype_usage, aes(fill = Isotype, y = Counts, x = Sample)) + geom_bar(position="stack", stat="identity") +ggtitle("Isotype usage of public cells")+ xlab("") + ylab("") + theme(plot.title = element_text(hjust = 0.5))

    #Look at scaled isotype usage (by number of public cells per sample)
    isotype_frequency_public_cells <- ggplot(isotype_usage, aes(fill = Isotype, y = Scaled_Counts, x = Sample)) + geom_bar(position="stack", stat="identity") +ggtitle("Isotype usage of public cells - scaled by # public cells within sample")+ xlab("") + ylab("") + theme(plot.title = element_text(hjust = 0.5))

    if(save.plot == TRUE){
      ggsave(isotype_usage_public_cells, filename = "isotype_usage_public_cells.png", dpi = 400, width = 12.0, height = 9.48)
      ggsave(isotype_frequency_public_cells, filename = "isotype_frequency_public_cells.png", dpi = 400, width = 12.0, height = 9.48)
    }

    tmp <- list(isotype_usage_public_cells, isotype_frequency_public_cells )
    names(tmp) <- c("isotype_usage_public_cells", "isotype_frequency_public_cells" )

    return(tmp)


  }


  # 6. Transcriptional cluster usage of public cells
  if(analysis.type == "transcriptional cluster usage of public cells"){

    #Create subset containing only cells with a public clonotype
    public_cells <- subset(vdj, vdj$public == TRUE)

    #replace all NA's with Unknown in the seurat cluster column
    public_cells$seurat_clusters <- factor(public_cells$seurat_clusters, levels = c(levels(public_cells$seurat_clusters), "Unknown"))
    public_cells$seurat_clusters[is.na(public_cells$seurat_clusters)] <- "Unknown"

    cluster_usage <- NULL

    for(i in 1:length(sample_ids)){
      curr_sample <- subset(public_cells, public_cells$sample_id == sample_ids[i])

      df <- as.data.frame(table(curr_sample$seurat_clusters))
      df$scaled_freq <- df$Freq/nrow(curr_sample)*100
      df$sample <- sample_ids[i]

      names(df) <- c("Cluster", "Counts", "Scaled_Counts", "Sample")

      cluster_usage <- rbind(cluster_usage,df)

    }


    cluster_usage_public_cells <- ggplot(cluster_usage, aes(fill = Cluster, y = Counts, x = Sample)) + geom_bar(position="stack", stat="identity") +ggtitle("Transcriptional clusters of public cells")+ xlab("") + ylab("") + theme(plot.title = element_text(hjust = 0.5))

    cluster_frequency_public_cells <- ggplot(cluster_usage, aes(fill = Cluster, y = Scaled_Counts, x = Sample)) + geom_bar(position="stack", stat="identity") +ggtitle("Transcriptional cluster frequency")+ xlab("") + ylab("") + theme(plot.title = element_text(hjust = 0.5))


    if(save.plot == TRUE){
      ggsave(cluster_usage_public_cells, filename = "cluster_usage_public_cells.png", dpi = 400, width = 12.0, height = 9.48)
      ggsave(cluster_frequency_public_cells, filename = "cluster_frequency_public_cells.png", dpi = 400, width = 12.0, height = 9.48)
    }

    tmp <- list(cluster_usage_public_cells, cluster_frequency_public_cells)
    names(tmp) <- c("cluster_usage_public_cells", "cluster_frequency_public_cells")

    return(tmp)


  }

}
