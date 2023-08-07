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

  platypus.version <- "v3"

  x <- NULL
  y <- NULL
  grob <- NULL
  width <- NULL
  height<- NULL

  ggname <- function(prefix, grob) {
    grob$name <- grid::grobName(grob, prefix)
    grob
  }

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
        ggname("geom_spatial", g)
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

  plot<-ggplot2::ggplot(data = vgm_VDJ, ggplot2::aes(x=x,y=y, fill = as.factor(analysis)))+
    geom_spatial(data=images_tibble[1,], ggplot2::aes(grob=grob), x=0.5, y=0.5)+
    ggplot2::geom_point(shape=21, colour = "black", size = 1.75, stroke = 0.5)+
    ggplot2::coord_cartesian(expand=FALSE)+
    ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(reverse=TRUE))+
    ggplot2::xlim(0,max(bcs_merge %>%
                          dplyr::filter(sample ==sample_names[1]) %>%
                          dplyr::select(width)))+
    ggplot2::ylim(max(bcs_merge %>%
                        dplyr::filter(sample ==sample_names[1]) %>%
                        dplyr::select(height)),0)+
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::ggtitle(sample_names[1], title)+
    ggplot2::theme(axis.text=ggplot2::element_text(size=size),
          axis.title=ggplot2::element_text(size=size))+
    ggplot2::labs(fill = legend_title)+
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size=3)))+
    ggplot2::theme_set(ggplot2::theme_bw(base_size = size))+
    ggplot2::theme(legend.key = ggplot2::element_rect(fill = "white"))+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank())
  return(plot)
}

