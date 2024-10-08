#' Plotting celltype assign to cell according to their phenotype on the spatial image.
#' @param sample_names Character vector containing the name of the sample.
#' @param bcs_merge Data frame containing imagerow, imagecol and barcode of the cells belonging to the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 10.
#' @param images_tibble Tbl-df containing the sample name, grob, height and width of the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 5.
#' @param vgm_GEX Data frame containing GEX information (VGM[[2]]). It must have a barcode column containing GEX_barcode and a cell.state column (output from GEX_phenotype).
#' @param title Character vector to name the plot.
#' @param size Number, to define the size of the text, default = 15.
#' @param legend_title Character vector to name the legend scale.
#' @param unclassified_cells Booleans, if TRUE the unclassified cells are also plot and if FALSE they aren't plot exept if the parameter specific_celltype = "Unclassified". In this case the unclassified cells are displayed even unclassified_cells = FALSE. Default = FALSE.
#' @param specific_celltype Character vector, the user can choose to express a specific celltype like T, B or Unclassified cells. Default = No.
#' @param density Booleans, if TRUE a density map is made. Default = FALSE
#' @return Returns a ggplot of the celltypes and if density = TRUE a density map of the cells on the spatial image.
#' @export
#' @examples
#' \dontrun{
#' Spatial_celltype_plot(bcs_merge = scaling_parameters[[10]],
#' vgm_GEX = vgm_spatial$GEX@meta.data,images_tibble = scaling_parameters[[5]],
#' sample_names = sample_names,title="B and T celltype", legend_title = "Celltype",
#' unclassified_cells = FALSE, specific_celltype = "Unclassified")
#' }

Spatial_celltype_plot<-function(sample_names,bcs_merge,images_tibble, vgm_GEX,  title, size, legend_title, unclassified_cells = c(TRUE, FALSE), specific_celltype = c("T","B","No", "Unclassified"),density=c(TRUE, FALSE)){
  if(missing(vgm_GEX)) stop("Please provide vgm_GEX input for this function")
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
    legend = ""
  }
  if (missing(unclassified_cells)){
    unclassified_cells = FALSE
  }
  if (missing(specific_celltype)){
    specific_celltype = "No"
  }
  if (missing(density)){
    density = FALSE
  }

  platypus.version <- "v3"

  x = NULL
  y = NULL
  cell.state = NULL
  grob = NULL
  celltype = NULL
  orig_barcode = NULL
  cell.state = NULL
  width = NULL
  height = NULL
  ..level.. = NULL

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

  GEX_celltype<-bcs_merge
  names(GEX_celltype)[6]<-"y"
  names(GEX_celltype)[7]<-"x"
  GEX_celltype$barcode<-gsub("-1","",as.character(GEX_celltype$barcode))
  celltype<-dplyr::select(vgm_GEX, orig_barcode, cell.state)
  names(celltype)[1]<-"barcode"
  GEX_celltype<-merge(GEX_celltype, celltype, by = "barcode")
  #specific T or B cell
  if(specific_celltype =="T"){
    GEX_celltype<-dplyr::filter(GEX_celltype, cell.state == "T")
  } else if (specific_celltype == "B"){
    GEX_celltype <-dplyr::filter(GEX_celltype, cell.state == "B")
  } else if (specific_celltype == "No"){
    GEX_celltype <- GEX_celltype
  } else if (specific_celltype =="Unclassified"){
    GEX_celltype <-dplyr::filter(GEX_celltype, cell.state == "Unclassified")
    unclassified_cells = TRUE
  }
  #with or without unclassified cells
  if (unclassified_cells == TRUE){
    GEX_celltype = GEX_celltype
  } else if(unclassified_cells == FALSE){
    GEX_celltype <- dplyr::filter(GEX_celltype, cell.state != "Unclassified")
  }

  plot<-ggplot2::ggplot(data = GEX_celltype, ggplot2::aes(x=x,y=y, fill = as.factor(cell.state)))+
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

  if (density == FALSE){
    return(plot)
  } else if (density == TRUE){
    density_plot<- ggplot2::ggplot(data = GEX_celltype, ggplot2::aes(x=x, y=y) ) +
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
    return(list(plot, density_plot))
  }
}



