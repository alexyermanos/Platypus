#' Plotting immune repertoire data as clonotype or isotype for cells on a spatial image.
#' @param sample_names Character vector containing the name of the sample.
#' @param bcs_merge Data frame containing imagerow, imagecol and barcode of the cells belonging to the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 10.
#' @param images_tibble Tbl-df containing the sample name, grob, height and width of the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 5.
#' @param title Character vector to name the plot.
#' @param size Number, to define the size of the text, default = 15.
#' @param legend_title Character vector to name the legend scale.
#' @param vgm_VDJ Data frame containing VDJ information, found in the vgm made by platypus. It must have x and y coordinates column and the column containing the factor to plot.
#' @param analysis Column in the dataframe containing the factor of interest to plot on the spatial image.
#' @return Returns a plot of the factor of interest express on a spatial image.
#' @export
#' @examples
#' \dontrun{
#' Spatial_VDJ_plot(vgm_VDJ = top_5_VDJ_data,analysis = top_5_VDJ_data$VDJ_cgene,
#' images_tibble = scaling_parameters[[5]], bcs_merge = scaling_parameters[[10]],
#' sample_names = sample_names, title = "B cell", legend = "Isotype")
#' }

Spatial_VDJ_plot<-function(sample_names,bcs_merge,images_tibble,title,size,legend_title,vgm_VDJ,analysis){
  
  if(missing(vgm_VDJ)) stop("Please provide vgm_VDJ input for this function")
  if(missing(analysis)) stop("Please provide analysis input for this function")
  if(missing(bcs_merge)) stop("Please provide bcs_merge input for this function")
  if(missing(images_tibble)) stop("Please provide images_tibble input for this function")
  if(missing(sample_names)) stop("Please provide sample_names input for this function")
  if(missing(size)){
    size = 15
  }
  if(missing(title)){
    title = ""
  }
  if(missing(legend_title)){
    legend_title = ""
  }
  ggplot(data = vgm_VDJ, aes(x=x,y=y, fill = as.factor(analysis)))+
    geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(shape=21, colour = "black", size = 1.75, stroke = 0.5)+
    coord_cartesian(expand=FALSE)+
    scale_fill_discrete(guide = guide_legend(reverse=TRUE))+
    xlim(0,max(bcs_merge %>% 
                 filter(sample ==sample_names[1]) %>% 
                 select(width)))+
    ylim(max(bcs_merge %>% 
               filter(sample ==sample_names[1]) %>% 
               select(height)),0)+
    xlab("") +
    ylab("") +
    ggtitle(sample_names[1], title)+
    theme(axis.text=element_text(size=size),
          axis.title=element_text(size=size))+
    labs(fill = legend_title)+
    guides(fill = guide_legend(override.aes = list(size=3)))+
    theme_set(theme_bw(base_size = size))+
    theme(legend.key = element_rect(fill = "white"))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

