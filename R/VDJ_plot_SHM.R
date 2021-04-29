#' Plots for SHM based on MIXCR output generated using the VDJ_call_MIXCR function and appended to the VDJ.GEX.matrix.output
#' @param VDJ.matrix VDJ dataframe generated using the VDJ_GEX_matrix function and supplemented with MIXCR information
#' @param group.by Character. Defaults to "sample_id". Column name of VDJ.matrix to split VDJ.matrix by. For each unique entry in that column a set of plots will be generated. This can be useful to plot SHM by expansion or by transcriptomics-derived clusters
#' @param quantile.label Numeric. Defaults to 0.9. Which points to label in the SHM scatterplot. If set to 0.9, the top 10% of cells by SHM number will be labelled. If ggrepel throws a warning, concerning overlap it is recommended to attempt to lable less points to avoid cluttering
#' @param platypus.version Character. Only "v3" available.  
#' @return Returns a list of ggplot objects. out[[1]] is a boxplot comparing SHM by group.by. out[[2]] to out[[n]] are plots for each group that visualize VDJ and VJ SHM distribution for each group. Data for any plot can be accessed via out[[any]]$data
#' @export
#' @examples
#' \dontrun{
#' 
#' Standard plot
#' SHM_plots <- plot_SHM(VDJ.matrix = VDJ_comb[[1]], group.by = "sample_id", quantile.label = 0.9)
#' 
#' Group by transcriptional cluster and label only top 1%
#' SHM_plots <- plot_SHM(VDJ.matrix = VDJ_comb[[1]], group.by = "seurat_clusters", quantile.label = 0.99)
#'}

plot_SHM <- function(VDJ.matrix, 
                     group.by, 
                     quantile.label,
                     platypus.version){
  
  
  
  platypus.version <- "v3"
  if(missing(quantile.label)) quantile.label <- 0.9
  if(missing(group.by)) group.by <- "sample_id"
  
  require(tidyr)
  require(ggrepel)
  
  #get data
  if(!"VDJ_SHM" %in% names(VDJ.matrix) | !"VJ_SHM" %in% names(VDJ.matrix)){
    stop("VDJ_SHM and VJ_SHM column not found in the input dataframe. Please provide a output dataframe of the VDJ_call_MIXCR function")
  }
  
  if((group.by %in% names(VDJ.matrix)) == F){
    stop(paste0("Specified group.by column " , group.by, " was not found in VDJ.matrix input dataframe. Please provide a valid column name"))
  }

  to_plot <- VDJ.matrix[,c(group.by, "barcode", "VDJ_SHM", "VJ_SHM")]
  names(to_plot)[1] <- "group"
  
  to_plot$VDJ_SHM[is.na(to_plot$VDJ_SHM)] <- 0
  to_plot$VJ_SHM[is.na(to_plot$VJ_SHM)] <- 0
  
  to_plot_long <- pivot_longer(to_plot, cols = c(3:4)) 
  
  #BOXPLOT COMPARING SHM PER GROUP
  
  box_plot <- ggplot2::ggplot(to_plot_long, ggplot2::aes(color = group, y= value, x= group)) + ggplot2::geom_boxplot(width=0.6) + geom_jitter(alpha = 0.5, width = 0.5) + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::ylab("SHM") + ggplot2::xlab("") + ggplot2::ggtitle(label = paste0("SHM per ", group.by))+ ggplot2::theme(strip.background = element_rect(color = "white", fill = "white")) + ggplot2::facet_wrap(~name) 
  #SCATTERPLOT FOR EVERY GROUP
  
  out.list <- list()
  out.list[[1]] <- box_plot
  
  for(j in 1:length(unique(to_plot$group))){
    
    curr_to_plot <- subset(to_plot, group == unique(to_plot$group)[j])

    qx_HC <- quantile(curr_to_plot$VDJ_SHM, probs = quantile.label)
    qx_LC <- quantile(curr_to_plot$VJ_SHM, probs = quantile.label)
    
    out.list[[j+1]] <- ggplot(curr_to_plot, aes(x = VDJ_SHM, y = VJ_SHM, col = VDJ_SHM + VJ_SHM)) + geom_point(show.legend = T, size = 3, alpha = 1) + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::ylab("VJ SHM") + ggplot2::xlab("VDJ SHM") + ggplot2::ggtitle(label = paste0("SHM in ", unique(to_plot$group)[j]))+ geom_text_repel(inherit.aes = T, data = subset(curr_to_plot, VDJ_SHM > qx_HC | VJ_SHM > qx_LC), aes(x = VDJ_SHM, y = VJ_SHM, label = barcode), color = "black") + scale_color_viridis_c(option = "B", end = 0.9)
    
  }
  
  return(out.list)
}
