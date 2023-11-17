#'Platypus V2 clonal utility
#'
#' @description Only for Platypus v2 Clonal frequency plot displaying the isotype usage of each clone. ! For platypus v3 use VDJ_clonal_expansion
#' @param VDJ_clonotype_output list of dataframes based on the VDJ_clonotype function output.
#' @param VDJ_per_clone_output list of dataframes based on the VDJ_per_clone function output.
#' @param clones numeric value indicating the number of clones to be displayed on the clonal expansion plot. Can take values between 1-50. Default value is 50.
#' @param subtypes Logical indicating whether to display isotype subtypes or not.
#' @param species Character indicating whether the samples are from mouse or human. Default is set to human.
#' #' @param sample.names Character vector with the same length of the VDJ.GEX.matrix.out list. If a VDJ table is provided, length of samples names must be one. These names are used as references to the output and as title for the plots
#' @param treat.incomplete.clones Character indicating how to proceed with clonotypes lacking a VDJC (in other words, no cell within the clonotype has a VDJC). "exclude" removes these clonotypes from the analysis. This may result in a different frequency ranking of clonotypes than in the output of the VDJ_analyse function with filter.1HC.1LC = FALSE. "include" keeps these clonotypes in the analysis. In the plot they will appear has having an unknown isotype.
#' @param treat.incomplete.cells Character indicating how to proceed with cells assigned to a clonotype but missing a VDJC. "proportional" to fill in the VDJ isotype according to the proportions present in of clonotype (in case present proportions are not replicable in the total number of cells e.g. 1/3 in 10 cells, values are rounded to the next full integer and if the new counts exceed the total number of cells, 1 is subtracted from the isotype of highest frequency. If the number is below the number of cell, 1 is added to the isotype with lowest frequency to preserve diversity), "exclude" to exclude them from analysis and rank clonotypes only by the number of actual contigs of there heavy chain. This ranking may deviate from the frequency column in the clonotype table. CAVE: if treat_incomplete_cells is set to "exclude", clonotypes lacking a VDJC entierly will be removed from the analysis. This results in a similar but not identical output as when treat_incomplete_clones is set to true. The two parameters are thereby non-redundant.
#' @param sample.names Vector. Names for samples in the order of the VDJ_GEX_matrix or the VDJ.analyze.output. Defaults to 1-n
#' @param platypus.version Defaults to "v3". For a more flexible analysis in v3 use VDJ_clonal_expansion()
#' @param VDJ.matrix The VDJ table output of the VDJ_GEX_matrix function. (VDJ_GEX_matrix.output[[1]])
#' @return returns a list containing plots with the percentages of isotypes for each clone on the cell level.
#' @export
#' @examples
#' \dontrun{
#' VDJ.isotype.per.clone <- VDJ_isotypes_per_clone(
#' VDJ_clonotype_output = VDJ.analyze.output
#' ,VDJ_per_clone_output = VDJ.per.clone.output, clones = 30)
#'}
VDJ_isotypes_per_clone <- function(VDJ_clonotype_output,
                                   VDJ_per_clone_output,
                                   clones,
                                   subtypes,
                                   species,
                                   sample.names,
                                   treat.incomplete.clones,
                                   treat.incomplete.cells,
                                   platypus.version,
                                   VDJ.matrix){

  Counts <- NULL
  sum_counts <- NULL
  Isotype <- NULL
  ClonalRank <- NULL


  if(missing(clones)) message("Number of clones to be displayed has not been supplied. 50 clones will be displayed by default")
  if(missing(clones)) clones <- 50
  if(missing(VDJ_clonotype_output)) VDJ_clonotype_output <- list()
  if(missing(VDJ_per_clone_output)) VDJ_per_clone_output <- list()
  if(missing(VDJ.matrix)) VDJ.matrix <- list()

  VDJ.GEX.matrix.out <- VDJ.matrix

  if(missing(subtypes)) subtypes <- FALSE
  if(missing(species)) species <- "Human"
  if(missing(treat.incomplete.cells)) treat.incomplete.cells <- "proportional"
  if(missing(treat.incomplete.clones)) treat.incomplete.clones <- "exclude"
  if(missing(sample.names)) sample.names <- c(1:length(VDJ.GEX.matrix.out))

  if(missing(platypus.version)) platypus.version <- "v3"

  if(platypus.version=="v3"){ ####START v3

      if(inherits(VDJ.GEX.matrix.out,"data.frame")){ #converting a single dataframe intput into a list for compatibility
        VDJ.GEX.matrix.out <- list(VDJ.GEX.matrix.out)
      }
    if(missing(sample.names)) sample.names <- c(1:length(VDJ.GEX.matrix.out))

    if(clones<1 | clones>50){stop("Number of clones must be an integer value between 1 and 50")}

    VDJ_per_clone_output_all <- list()
    clones_per_isotype <- list()
    clones_per_isotype_all <-list()
    output_plot <- list()

    if(subtypes == FALSE){

      for (i in 1:length(VDJ.GEX.matrix.out)){
        #get clonotype frequency table
        clono_freq <- as.data.frame(table(VDJ.GEX.matrix.out[[i]]$clonotype_id_10x))
        clono_freq <- clono_freq[order(clono_freq$Freq, decreasing = T),]

        #get essential info from VDJ_GEX_matrix
        curr_rep_iso <- VDJ.GEX.matrix.out[[i]][,c("barcode","clonotype_id_10x", "VDJ_cgene", "VDJ_cdr3s_aa", "VJ_cdr3s_aa")]
        curr_rep_iso$isotype <- substring(curr_rep_iso$VDJ_cgene,first = 1,last = 4)

        if(treat.incomplete.clones == "exclude"){
          clones_to_del <-c()
          for(k in 1:nrow(clono_freq)){
            if(all(curr_rep_iso[which(curr_rep_iso$clonotype_id_10x == clono_freq[k,1]),"VDJ_cgene"] == "")){ #detect if no VDJC is available in this clonotype
              clones_to_del <- append(clones_to_del, k)
            }
          }
          if(length(clones_to_del > 0)){
            clono_freq <- clono_freq[-clones_to_del,]} #remove these clones entirely
          clono_freq <- clono_freq[order(clono_freq$Freq, decreasing = T),] #reorder, to make sure
        }

        #iterate over clones and get isotype info
        clones_per_isotype <- list()
        for (j in 1:clones){
          curr_clone <- curr_rep_iso[which(curr_rep_iso$clonotype_id_10x == clono_freq[j,1]),]

          if(treat.incomplete.cells == "proportional" & stringr::str_detect(paste0(curr_clone$isotype,collapse = ";"), pattern = "IG")==T){ #check that there is at least one IG entry
            props <- table(curr_clone$isotype[which(stringr::str_detect(curr_clone$isotype, pattern = "IG")==T)]) #getting proportions
            n_total <- nrow(curr_clone) #getting total number of cells in clonotype
            props <- round(props / sum(props) * n_total,0) #calculating new number of each isotype to match proportions
            if(sum(props) > n_total){props[which.max(props)] <- props[which.max(props)] - 1
            } else if(sum(props) < n_total){props[which.min(props)] <- props[which.min(props)] + 1} #catching rounding derived errors.

            curr_clone$isotype <- rep.int(names(props), props) #new isotype column with the new number of isotypes. ! these are not in the original order !

          } else if (treat.incomplete.cells == "proportional" & stringr::str_detect(paste0(curr_clone$isotype,collapse = ";"), pattern = "IG")==F){ #if no entry is present

            curr_clone$isotype <- "None"
          }

          clones_per_isotype[[j]] <- data.frame("Counts"=rep(0, 6), "Color"=rep("", 6), "Isotype"=rep("", 6), "ClonalRank"=rep("", 6), "clonotype_id" = rep(curr_clone$clonotype_id_10x[1],6), "VDJ_cdr3s_aa" = rep(curr_clone$VDJ_cdr3s_aa[which(curr_clone$VDJ_cdr3s_aa != "")][1],6), "VJ_cdr3s_aa" = rep(curr_clone$VJ_cdr3s_aa[which(curr_clone$VJ_cdr3s_aa != "")][1],6), "barcode" = rep(paste0(curr_clone$barcode, collapse = ";"),6))  #to maintain clonotype information

          clones_per_isotype[[j]]$Counts[1] <- sum(stringr::str_count(curr_clone$isotype, "IGHG"))
          clones_per_isotype[[j]]$Counts[2] <- sum(stringr::str_count(curr_clone$isotype, "IGHM"))
          clones_per_isotype[[j]]$Counts[3] <- sum(stringr::str_count(curr_clone$isotype, "IGHA"))
          clones_per_isotype[[j]]$Counts[4] <- sum(stringr::str_count(curr_clone$isotype, "IGHD"))
          clones_per_isotype[[j]]$Counts[5] <- sum(stringr::str_count(curr_clone$isotype, "IGHE"))
          clones_per_isotype[[j]]$Counts[6] <- sum(stringr::str_count(curr_clone$isotype, "None"))

          clones_per_isotype[[j]]$Color <- c("green4", "black", "red", "blue", "purple", "gray")
          clones_per_isotype[[j]]$Isotype <- c("IGHG", "IGHM", "IGHA", "IGHD", "IGHE", "Unknown")
          clones_per_isotype[[j]]$ClonalRank <- j


        }
        clones_per_isotype_all[[i]] <- do.call("rbind", clones_per_isotype)

        if(treat.incomplete.cells == "exclude"){
          rank_raw <- as.data.frame(clones_per_isotype_all[[i]] %>% dplyr::group_by(ClonalRank) %>% dplyr::summarise(sum_counts = sum(Counts)) %>% dplyr::arrange(dplyr::desc(sum_counts)) %>% dplyr::mutate(rank = 1:length(unique(ClonalRank))))
          clones_per_isotype_all[[i]]$ClonalRank_2 <- 0
          for(l in 1:nrow(rank_raw)){
            clones_per_isotype_all[[i]]$ClonalRank_2[which(clones_per_isotype_all[[i]]$ClonalRank == rank_raw$ClonalRank[l])] <- rank_raw$rank[l]
          }
          clones_per_isotype_all[[i]]$ClonalRank <- clones_per_isotype_all[[i]]$ClonalRank_2
          #reorder the dataframe the sum_counts order by decreasing sum, the clonal rank is there to keep clones that have the same count sum together as a group of rows

          message("New ranking based only on present HC chains: ")
          message(unique(clones_per_isotype_all[[i]]$clonotype_id))
        }

        output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], ggplot2::aes(fill = Isotype, y=Counts, x=ClonalRank)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black") + ggplot2::theme_bw() + ggplot2::scale_fill_manual("Isotype", values = c("IGHG" = "green4", "IGHM" = "black", "IGHA" = "red3", "IGHD"="blue", "IGHE"="purple", "Unknown"="gray")) + ggplot2::theme_classic() + ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::scale_x_continuous(expand = c(0,0.5)) + ggplot2::labs(title = sample.names[[i]], x = "Clonal rank", y = "Number of cells")

      }
      names(clones_per_isotype_all) <- sample.names
      return(list(clones_per_isotype_all,output_plot))
    }
    if(subtypes == TRUE){

      for (i in 1:length(VDJ.GEX.matrix.out)){
        #get clonotype frequency table
        clono_freq <- as.data.frame(table(VDJ.GEX.matrix.out[[i]]$clonotype_id_10x))
        clono_freq <- clono_freq[order(clono_freq$Freq, decreasing = T),]

        #get essential info from VDJ_GEX_matrix
        curr_rep_iso <- VDJ.GEX.matrix.out[[i]][,c("barcode","clonotype_id_10x", "VDJ_cgene", "VDJ_cdr3s_aa", "VJ_cdr3s_aa")]
        #no substring here, just splitting on a ; to exclude multiple isotypes of one clone
        curr_rep_iso$isotype <- stringr::str_split(curr_rep_iso$VDJ_cgene, ";", simplify = T)[,1]



        if(treat.incomplete.clones == "exclude"){
          clones_to_del <-c()
          for(k in 1:nrow(clono_freq)){
            if(all(curr_rep_iso[which(curr_rep_iso$clonotype_id_10x == clono_freq[k,1]),"VDJ_cgene"] == "")){ #detect if no VDJC is available in this clonotype
              clones_to_del <- append(clones_to_del, k)
            }
          }
          if(length(clones_to_del > 0)){
            clono_freq <- clono_freq[-clones_to_del,]} #remove these clones entirely
          clono_freq <- clono_freq[order(clono_freq$Freq, decreasing = T),] #reorder, to make sure
        }

        #iterate over clones and get isotype info
        clones_per_isotype <- list()
        j <- 1
        for (j in 1:clones){
          curr_clone <- curr_rep_iso[which(curr_rep_iso$clonotype_id_10x == clono_freq[j,1]),]

          if(treat.incomplete.cells == "proportional" & stringr::str_detect(paste0(curr_clone$isotype,collapse = ";"), pattern = "IG")==T){ #check that there is at least one IG entry
            props <- table(curr_clone$isotype[which(stringr::str_detect(curr_clone$isotype, pattern = "IG")==T)]) #getting proportions
            n_total <- nrow(curr_clone) #getting total number of cells in clonotype
            props <- round(props / sum(props) * n_total,0) #calculating new number of each isotype to match proportions
            if(sum(props) > n_total){props[which.max(props)] <- props[which.max(props)] - 1
            } else if(sum(props) < n_total){props[which.min(props)] <- props[which.min(props)] + 1} #catching rounding derived errors.

            curr_clone$isotype <- rep.int(names(props), props) #new isotype column with the new number of isotypes. ! these are not in the original order !

          } else if (treat.incomplete.cells == "proportional" & stringr::str_detect(paste0(curr_clone$isotype,collapse = ";"), pattern = "IG")==F){ #if no entry is present

            curr_clone$isotype <- "None"
          }

          clones_per_isotype[[j]] <- data.frame("Counts"=rep(0, 14), "Color"=rep("", 14), "Isotype"=rep("", 14), "ClonalRank"=rep("", 14), "clonotype_id" = rep(curr_clone$clonotype_id_10x[1],14), "VDJ_cdr3s_aa" = rep(curr_clone$VDJ_cdr3s_aa[which(curr_clone$VDJ_cdr3s_aa != "")][1],14), "VJ_cdr3s_aa" = rep(curr_clone$VJ_cdr3s_aa[which(curr_clone$VJ_cdr3s_aa != "")][1],14), "barcode" = rep(paste0(curr_clone$barcode, collapse = ";"),14))  #to maintain clonotype information

          clones_per_isotype[[j]]$Counts[1] <- sum(stringr::str_count(curr_clone$isotype, "IGHG1"))
          clones_per_isotype[[j]]$Counts[2] <- sum(stringr::str_count(curr_clone$isotype, "IGHG2"))
          clones_per_isotype[[j]]$Counts[3] <- sum(stringr::str_count(curr_clone$isotype, "IGHG2A"))
          clones_per_isotype[[j]]$Counts[4] <- sum(stringr::str_count(curr_clone$isotype, "IGHG2B"))
          clones_per_isotype[[j]]$Counts[5] <- sum(stringr::str_count(curr_clone$isotype, "IGHG2C"))
          clones_per_isotype[[j]]$Counts[6] <- sum(stringr::str_count(curr_clone$isotype, "IGHG3"))
          clones_per_isotype[[j]]$Counts[7] <- sum(stringr::str_count(curr_clone$isotype, "IGHG4"))
          clones_per_isotype[[j]]$Counts[8] <- sum(stringr::str_count(curr_clone$isotype, "IGHM"))
          clones_per_isotype[[j]]$Counts[9] <- sum(stringr::str_count(curr_clone$isotype, "IGHA"))
          clones_per_isotype[[j]]$Counts[10] <- sum(stringr::str_count(curr_clone$isotype, "IGHA1"))
          clones_per_isotype[[j]]$Counts[11] <- sum(stringr::str_count(curr_clone$isotype, "IGHA2"))
          clones_per_isotype[[j]]$Counts[12] <- sum(stringr::str_count(curr_clone$isotype, "IGHD"))
          clones_per_isotype[[j]]$Counts[13] <- sum(stringr::str_count(curr_clone$isotype, "IGHE"))
          clones_per_isotype[[j]]$Counts[14] <- sum(stringr::str_count(curr_clone$isotype, "None"))


          clones_per_isotype[[j]]$Isotype <- c("IGHG1", "IGHG2", "IGHG2a", "IGHG2b", "IGHG2c", "IGHG3", "IGHG4", "IGHM", "IGHA", "IGHA1", "IGHA2", "IGHD", "IGHE", "Unknown")
          clones_per_isotype[[j]]$ClonalRank <- j

          if (species == "Human"){
            clones_per_isotype[[j]] <- clones_per_isotype[[j]][-c(3, 4, 5, 9), ]
          }

          if (species == "Mouse"){
            clones_per_isotype[[j]] <- clones_per_isotype[[j]][-c(2, 7, 10, 11), ]
          }


        }
        clones_per_isotype_all[[i]] <- do.call("rbind",clones_per_isotype)

        if(treat.incomplete.cells == "exclude"){
          rank_raw <- as.data.frame(clones_per_isotype_all[[i]] %>% dplyr::group_by(ClonalRank) %>% dplyr::summarise(sum_counts = sum(Counts)) %>% dplyr::arrange(dplyr::desc(sum_counts)) %>% dplyr::mutate(rank = 1:length(unique(ClonalRank))))
          clones_per_isotype_all[[i]]$ClonalRank_2 <- 0
          for(l in 1:nrow(rank_raw)){
            clones_per_isotype_all[[i]]$ClonalRank_2[which(clones_per_isotype_all[[i]]$ClonalRank == rank_raw$ClonalRank[l])] <- rank_raw$rank[l]
          }

          clones_per_isotype_all[[i]]$ClonalRank <- clones_per_isotype_all[[i]]$ClonalRank_2
          #reorder the dataframe the sum_counts order by decreasing sum, the clonal rank is there to keep clones that have the same count sum together as a group of rows

          message("New ranking based only on present HC chains: ")
          message(unique(clones_per_isotype_all[[i]]$clonotype_id))

        }
        if (species == "Human"){
          output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], ggplot2::aes(fill = Isotype, y=Counts, x=ClonalRank)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black") + cowplot::theme_cowplot() + ggplot2::scale_fill_manual("Isotype", values = c("IGHG1" = "green","IGHG2" = "green3", "IGHG3"="green4", "IGHG4"="darkgreen", "IGHM" = "black", "IGHA1" = "red", "IGHA2"= "red4", "IGHD"="blue", "IGHE"="purple", "Unknown"="gray")) + ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Clonal rank") + ggplot2::scale_x_continuous(expand = c(0,0.5))
        }

        if (species == "Mouse"){
          output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], ggplot2::aes(fill = Isotype, y=Counts, x=ClonalRank)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black") + cowplot::theme_cowplot() + ggplot2::scale_fill_manual("Isotype", values = c("IGHG1" = "lightgreen", "IGHG2a"="green", "IGHG2b" = "green3", "IGHG2c"="green4", "IGHG3"="darkgreen", "IGHM" = "black", "IGHA" = "red", "IGHD"="blue", "IGHE"="purple", "Unknown"="gray"))  + ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Clonal rank") + ggplot2::scale_x_continuous(expand = c(0,0.5))

        }
      }
    }

    return(list(clones_per_isotype_all,output_plot))
  }####STOP v3
  else if(platypus.version=="v2"){ ####START v2



  if(clones<1 | clones>50){
    stop("Number of clones must be an integer value between 1 and 50")
  }else{

    VDJ_per_clone_output_all <- list()
    clones_per_isotype <- list()
    clones_per_isotype_all <-list()
    output_plot <- list()

    if(subtypes == FALSE){

      for (i in 1:length(VDJ_per_clone_output)){

        for (j in 1:length(VDJ_per_clone_output[[i]])){
          VDJ_per_clone_output[[i]][[j]] <- VDJ_per_clone_output[[i]][[j]][,c("barcode", "isotype_hc")]
        }

        VDJ_per_clone_output_all[[i]] <- do.call("rbind",VDJ_per_clone_output[[i]])
        VDJ_per_clone_output_all[[i]]$isotype <- substring(VDJ_per_clone_output_all[[i]]$isotype_hc,first = 1,last = 4)

        #Order the VDJ_clonotype output based on increasing frequency
        ranks <- order(-VDJ_clonotype_output[[i]]$frequency)
        VDJ_clonotype_output[[i]] <- VDJ_clonotype_output[[i]][ranks,]

        for (j in 1:clones){
          VDJ_clonotype_output_split <- as.data.frame(stringr::str_split(string = VDJ_clonotype_output[[i]]$barcodes[j],pattern = ";")[[1]])
          colnames(VDJ_clonotype_output_split) <- "barcode"

          for (k in 1:nrow(VDJ_clonotype_output_split)){
            VDJ_clonotype_output_split$isotype[k] <- VDJ_per_clone_output_all[[i]]$isotype[which(VDJ_per_clone_output_all[[i]]$barcode == VDJ_clonotype_output_split$barcode[k])]
          }

          clones_per_isotype[[j]] <- data.frame("Counts"=rep(0, 6), "Color"=rep("", 6), "Isotype"=rep("", 6), "ClonalRank"=rep("", 6))

          clones_per_isotype[[j]]$Counts[1] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHG"))
          clones_per_isotype[[j]]$Counts[2] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHM"))
          clones_per_isotype[[j]]$Counts[3] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHA"))
          clones_per_isotype[[j]]$Counts[4] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHD"))
          clones_per_isotype[[j]]$Counts[5] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHE"))
          clones_per_isotype[[j]]$Counts[6] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "None"))

          clones_per_isotype[[j]]$Color <- c("green4", "black", "red", "blue", "purple", "gray")
          clones_per_isotype[[j]]$Isotype <- c("IGHG", "IGHM", "IGHA", "IGHD", "IGHE", "Unknown")
          clones_per_isotype[[j]]$ClonalRank <- j

        }
        clones_per_isotype_all[[i]] <- do.call("rbind",clones_per_isotype)
        output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], ggplot2::aes(fill = Isotype, y=Counts, x=ClonalRank)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black") + cowplot::theme_cowplot() + ggplot2::scale_fill_manual("Isotype", values = c("IGHG" = "green4", "IGHM" = "black", "IGHA" = "red3", "IGHD"="blue", "IGHE"="purple", "Unknown"="gray")) + ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Clonal rank") + ggplot2::scale_x_continuous(expand = c(0,0.5))

      }
    }
      if(subtypes == TRUE){

        for (i in 1:length(VDJ_per_clone_output)){

          for (j in 1:length(VDJ_per_clone_output[[i]])){
            VDJ_per_clone_output[[i]][[j]] <- VDJ_per_clone_output[[i]][[j]][,c("barcode", "isotype_hc")]
          }

          VDJ_per_clone_output_all[[i]] <- do.call("rbind",VDJ_per_clone_output[[i]])

          #Order the VDJ_clonotype output based on increasing frequency
          ranks <- order(-VDJ_clonotype_output[[i]]$frequency)
          VDJ_clonotype_output[[i]] <- VDJ_clonotype_output[[i]][ranks,]

          for (j in 1:clones){
            VDJ_clonotype_output_split <- as.data.frame(stringr::str_split(string = VDJ_clonotype_output[[i]]$barcodes[j],pattern = ";")[[1]])
            colnames(VDJ_clonotype_output_split) <- "barcode"

            for (k in 1:nrow(VDJ_clonotype_output_split)){
              VDJ_clonotype_output_split$isotype[k] <- VDJ_per_clone_output_all[[i]]$isotype[which(VDJ_per_clone_output_all[[i]]$barcode == VDJ_clonotype_output_split$barcode[k])]
            }

            clones_per_isotype[[j]] <- data.frame("Counts"=rep(0, 14), "Isotype"=rep("", 14), "ClonalRank"=rep("", 14))

            clones_per_isotype[[j]]$Counts[1] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHG1"))
            clones_per_isotype[[j]]$Counts[2] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHG2"))
            clones_per_isotype[[j]]$Counts[3] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHG2a"))
            clones_per_isotype[[j]]$Counts[4] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHG2b"))
            clones_per_isotype[[j]]$Counts[5] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHG2c"))
            clones_per_isotype[[j]]$Counts[6] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHG3"))
            clones_per_isotype[[j]]$Counts[7] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHG4"))
            clones_per_isotype[[j]]$Counts[8] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHM"))
            clones_per_isotype[[j]]$Counts[9] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHA"))
            clones_per_isotype[[j]]$Counts[10] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHA1"))
            clones_per_isotype[[j]]$Counts[11] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHA2"))
            clones_per_isotype[[j]]$Counts[12] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHD"))
            clones_per_isotype[[j]]$Counts[13] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "IGHE"))
            clones_per_isotype[[j]]$Counts[14] <- sum(stringr::str_count(VDJ_clonotype_output_split$isotype, "None"))


            clones_per_isotype[[j]]$Isotype <- c("IGHG1", "IGHG2", "IGHG2a", "IGHG2b", "IGHG2c", "IGHG3", "IGHG4", "IGHM", "IGHA", "IGHA1", "IGHA2", "IGHD", "IGHE", "Unknown")
            clones_per_isotype[[j]]$ClonalRank <- j

            if (species == "Human"){
              clones_per_isotype[[j]] <- clones_per_isotype[[j]][-c(3, 4, 5, 9), ]
            }

            if (species == "Mouse"){
              clones_per_isotype[[j]] <- clones_per_isotype[[j]][-c(2, 7, 10, 11), ]
            }
          }

          clones_per_isotype_all[[i]] <- do.call("rbind",clones_per_isotype)

          if (species == "Human"){
          output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], ggplot2::aes(fill = Isotype, y=Counts, x=ClonalRank)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black") + cowplot::theme_cowplot() + ggplot2::scale_fill_manual("Isotype", values = c("IGHG1" = "green","IGHG2" = "green3", "IGHG3"="green4", "IGHG4"="darkgreen", "IGHM" = "black", "IGHA1" = "red", "IGHA2"= "red4", "IGHD"="blue", "IGHE"="purple", "Unknown"="gray")) + ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Clonal rank") + ggplot2::scale_x_continuous(expand = c(0,0.5))
          }

          if (species == "Mouse"){
          output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], ggplot2::aes(fill = Isotype, y=Counts, x=ClonalRank)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black") + cowplot::theme_cowplot() + ggplot2::scale_fill_manual("Isotype", values = c("IGHG1" = "lightgreen", "IGHG2a"="green", "IGHG2b" = "green3", "IGHG2c"="green4", "IGHG3"="darkgreen", "IGHM" = "black", "IGHA" = "red", "IGHD"="blue", "IGHE"="purple", "Unknown"="gray")) + ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Clonal rank") + ggplot2::scale_x_continuous(expand = c(0,0.5))
          }
        }
      }

    return(output_plot)
  }
  }####STOP v2
}
