#' Plotting the contour density of selected cells or of all cells.
#' @param sample_names Character vector containing the name of the sample.
#' @param bcs_merge Data frame containing imagerow, imagecol and barcode of the cells belonging to the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 10.
#' @param images_tibble Tbl-df containing the sample name, grob, height and width of the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 5.
#' @param vgm_VDJ Data frame containing all the data on the cell. It must contain the column clonotype_id which describes the number of the clonotype to which the cell belongs. This data frame can be obtained by the assignment functions (VDJ_assignment_random_based, VDJ_assignment_density_based and VDJ_assignment_germline_based).
#' @param title Character vector to name the plot.
#' @param size Number, to define the size of the text, default = 15.
#' @return Returns a plot of cell contour density on the spatial image.
#' @export
#' @examples
#' \dontrun{
#' #Assignment density-based
#' density_BCR_assignment<-Spatial_VDJ_assignment(GEX_matrix = GEX_matrix,vgm = vgm_with_simulated_VDJ,
#' vgm_VDJ = vgm_with_simulated_VDJ$VDJ, celltype = "B", simulated_VDJ = simulated_B_cells_VDJ, 
#' method = "density")
#' vgm_with_simulated_VDJ$VDJ<-density_BCR_assignment
#' 
#' top_1_VDJ_BCR_density_data<-Spatial_selection_expanded_clonotypes(nb_clonotype = 1, vgm_VDJ = vgm_with_simulated_VDJ$VDJ)
#' p_spatial_BCR_density_clonotype_density<-Spatial_density_plot(vgm_VDJ = top_1_VDJ_BCR_density_data,
#' images_tibble = scaling_parameters[[5]], bcs_merge = scaling_parameters[[10]],sample_names = sample_names,
#' title = "B cell density assignment")
#' p_spatial_BCR_density_clonotype_density
#' }
Spatial_density_plot<-function(sample_names,bcs_merge,images_tibble,vgm_VDJ,title,size){
  
  if(missing(sample_names)) stop("Please provide sample_names input for this function")
  if(missing(size)){
    size = 15
  }
  if(missing(images_tibble)) stop("Please provide images_tibble input for this function")
  if(missing(bcs_merge)) stop("Please provide bcs_merge input for this function")
  if(missing(vgm_VDJ)) stop("Please provide vgm_VDJ input for this function")
  if(missing(title)){
    title <- ""}
  
  ggplot(data = vgm_VDJ, aes(x=x, y=y) ) +
    geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)+
    coord_cartesian(expand=FALSE)+
    stat_density_2d(aes(fill = ..level..), alpha = 0.2, geom = "polygon", colour="white")+
    scale_fill_viridis_c()+
    xlim(0,max(bcs_merge %>% 
                 filter(sample ==sample_names[1]) %>% 
                 select(width)))+
    ylim(max(bcs_merge %>% 
               filter(sample ==sample_names[1]) %>% 
               select(height)),0)+
    xlab("") +
    ylab("") +
    ggtitle(sample_names[1],title)+
    theme(axis.text=element_text(size=size),
          axis.title=element_text(size=size))+
    labs(fill = "Density")+
    theme_set(theme_bw(base_size = size))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

