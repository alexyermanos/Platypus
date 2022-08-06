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
#' density_BCR_assignment<-Spatial_VDJ_assignment(GEX_matrix = GEX_matrix,
#' vgm = vgm_with_simulated_VDJ,
#' vgm_VDJ = vgm_with_simulated_VDJ$VDJ, celltype = "B",
#' simulated_VDJ = simulated_B_cells_VDJ,
#' method = "density")
#' vgm_with_simulated_VDJ$VDJ<-density_BCR_assignment
#'
#' top_1_VDJ_BCR_density_data<-Spatial_selection_expanded_clonotypes(
#' nb_clonotype = 1, vgm_VDJ = vgm_with_simulated_VDJ$VDJ)
#'
#' p_spatial_BCR_density_clonotype_density<-Spatial_density_plot(
#' vgm_VDJ = top_1_VDJ_BCR_density_data,
#' images_tibble = scaling_parameters[[5]],
#' bcs_merge = scaling_parameters[[10]],sample_names = sample_names,
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

  platypus.version <- "v3"

  x = NULL
  y = NULL
  grob = NULL
  ..level..= NULL
  width = NULL
  height = NULL

  geom_spatial <-  function(mapping = NULL,
                            data = NULL,
                            stat = "identity",
                            position = "identity",
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = FALSE,
                            ...) {

    GeomCustom <- ggplot2::ggproto(
      "GeomCustom",
      ggplot2::Geom,
      setup_data = function(self, data, params) {
        data <- ggplot2::ggproto_parent(ggplot2::Geom, self)$setup_data(data, params)
        data
      },

      draw_group = function(data, panel_scales, coord) {
        vp <- grid::viewport(x=data$x, y=data$y)
        g <- grid::editGrob(data$grob[[1]], vp=vp)
        ggplot2:::ggname("geom_spatial", g)
      },

      required_aes = c("grob","x","y")

    )

    ggplot2::layer(
      geom = GeomCustom,
      mapping = mapping,
      data = data,
      stat = stat,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(na.rm = na.rm, ...)
    )
  }

  plot<-ggplot2::ggplot(data = vgm_VDJ, ggplot2::aes(x=x, y=y) ) +
    geom_spatial(data=images_tibble[1,], ggplot2::aes(grob=grob), x=0.5, y=0.5)+
    ggplot2::coord_cartesian(expand=FALSE)+
    ggplot2::stat_density_2d(ggplot2::aes(fill = ..level..), alpha = 0.2, geom = "polygon", colour="white")+
    ggplot2::scale_fill_viridis_c()+
    ggplot2::xlim(0,max(bcs_merge %>%
                 dplyr::filter(sample ==sample_names[1]) %>%
                   dplyr::select(width)))+
    ggplot2::ylim(max(bcs_merge %>%
                        dplyr::filter(sample ==sample_names[1]) %>%
                        dplyr::select(height)),0)+
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::ggtitle(sample_names[1],title)+
    ggplot2::theme(axis.text=ggplot2::element_text(size=size),
          axis.title=ggplot2::element_text(size=size))+
    ggplot2::labs(fill = "Density")+
    ggplot2::theme_set(ggplot2::theme_bw(base_size = size))+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank())
  return(plot)
}

