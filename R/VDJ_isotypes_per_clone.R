#' Clonal frequency plot displaying the isotype usage of each clone.
#' @param VDJ_clonotype_output list of dataframes based on the VDJ_clonotype function output.
#' @param VDJ_per_clone_output list of dataframes based on the VDJ_per_clone function output.
#' @param clones numeric value indicating the number of clones to be displayed on the clonal expansion plot. Can take values between 1-50. Default value is 50.
#' @return returns a list containing plots with the percentages of isotypes for each clone on the cell level.
#' @export
#' @examples
#' \dontrun{
#' VDJ.isotype.per.clone <- VDJ_isotypes_per_clone(VDJ_clonotype_output = VDJ.analyze.output, VDJ_per_clone_output = VDJ.per.clone.output, clones = 30)
#'}
VDJ_isotypes_per_clone <- function(VDJ_clonotype_output,
                                   VDJ_per_clone_output,
                                   clones){
  require(stringr)
  require(ggplot2)

  if(missing(clones)) print("Number of clones to be displayed has not been supplied. 50 clones will be displayed by default")
  if(missing(clones)) clones <- 50

  if(clones<1 | clones>50){
    print("Number of clones must be a numerical value between 1 and 50")
  }else{

    VDJ_per_clone_output_all <- list()
    clones_per_isotype <- list()
    clones_per_isotype_all <-list()
    output_plot <- list()

    for (i in 1:length(VDJ_per_clone_output)){

      for (j in 1:length(VDJ_per_clone_output[[i]])){
        VDJ_per_clone_output[[i]][[j]] <- VDJ_per_clone_output[[i]][[j]][,c("barcode", "c_gene_HC")]
      }

      VDJ_per_clone_output_all[[i]] <- do.call("rbind", VDJ_per_clone_output[[i]])
      VDJ_per_clone_output_all[[i]]$isotype <- substring(VDJ_per_clone_output_all[[i]]$c_gene_HC,first = 1,last = 4)

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

        clones_per_isotype[[j]]$Color <- c("green", "black", "red", "blue", "purple", "gray")
        clones_per_isotype[[j]]$Isotype <- c("IGHG", "IGHM", "IGHA", "IGHD", "IGHE", "None")
        clones_per_isotype[[j]]$ClonalRank <- j

      }
      clones_per_isotype_all[[i]] <- do.call("rbind",clones_per_isotype)
      output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], ggplot2::aes(fill = Isotype, y=Counts, x=ClonalRank)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black") + ggplot2::theme_bw() + ggplot2::scale_fill_manual("Isotype", values = c("IGHG" = "green", "IGHM" = "black", "IGHA" = "red", "IGHD"="blue", "IGHE"="purple", "None"="gray")) + ggplot2::theme_classic() + ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Clonal rank") + ggplot2::scale_x_continuous(expand = c(0,0.5))

    }

    return(output_plot)
  }
}
