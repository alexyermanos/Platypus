#' Integrates VDJ and gene expression libraries by providing cluster membership seq_per_vdj object. Output will plot which transcriptional cluster (GEX) that the cells of a given clonotype are found in.
#' @param GEX.list The output of the automate_GEX function.
#' @param VDJ.GEX.integrate.list Output from VDJ_GEX_integrate function. This object needs to have the GEX and VDJ information combined and integrated. This should be on the CLONAL level from the VDJ_GEX_integrate function.
#' @param highlight.isotype (Optional) isotype to plot, choose between ["None","A","E","M","G","G1","G2A","G2B","G2C","G3"]. Default is None.
#' @param highlight.number A vector corresponding to the rank of the clones that should be specified. Default is set to "20", which will present the cluster distribution for the top 20 clones.
#' @return ggplot2 plot that breaks down clonotype membership per cluster for the specified input clones.
#' @export
#' @examples
#' \dontrun{
#' vdj.gex.expansion <- VDJ_GEX_expansion(GEX.list=GEX.list.output[[1]]
#' ,VDJ.GEX.integrate.list=vdj.gex.integrate.output
#' ,highlight.isotype = "None",highlight.number=1:20)
#' }
VDJ_GEX_expansion <- function(GEX.list,
                              VDJ.GEX.integrate.list,
                              highlight.isotype,
                              highlight.number){

  cluster <- NULL
  value <- NULL
  L1 <- NULL

  if(missing(GEX.list)) stop("No provided list of Seurat objects. Please provide the output of the automate_GEX function.")
  if(missing(VDJ.GEX.integrate.list)) stop("No provided per clone VDJ breakdown. Please provide the output of the VDJ_per_clone function.")
  if(missing(highlight.isotype)) highlight.isotype <- "None"
  if(missing(highlight.number)) highlight.number <- 20

  if(!(mode(highlight.number)=="numeric")){
    stop("highlight.number can either be an integer or a vector of integers")
  }

  # initialize plot list
  pl_list <- list()
  for(i in 1:length(VDJ.GEX.integrate.list)){

    holding_clusters <- list()
    holding_all_clusters_integer <- list()
    for(j in 1:length(highlight.number)){
      if(highlight.isotype=="None"){
        holding_clusters[[j]] <- VDJ.GEX.integrate.list[[i]]$cluster_membership_percent[highlight.number[j]]
        holding_all_clusters_integer[[j]] <- as.numeric(stringr::str_split(holding_clusters[[j]],pattern = ", ")[[1]])
        names(holding_all_clusters_integer[[j]]) <- paste("cluster",0:(length(holding_all_clusters_integer[[j]])-1),sep="")
      }
    }
    temp_melt <- reshape2::melt(holding_all_clusters_integer)
    temp_melt$cluster <- names(unlist(holding_all_clusters_integer))
    temp_melt$cluster <-  factor(temp_melt$cluster, levels = (names(holding_all_clusters_integer[[1]])))
    pl_list[[i]] <- ggplot2::ggplot(temp_melt, ggplot2::aes(fill = cluster, y=value, x=as.integer(L1))) + ggplot2::geom_bar(stat="identity", width=0.6, color="black") + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Percentage of cells") + ggplot2::xlab("Clonal rank") + ggplot2::scale_x_continuous(expand = c(0,0.5), breaks = c(1:length(highlight.number)))
  }
  return(pl_list)
}
