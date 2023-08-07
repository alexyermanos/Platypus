#' Plotting clusters of cells by choosing between 10X Genomics clustering or reclustering the cells.
#' @param cluster Character vector to describe the clustering, "GEX_cluster" is for plotting 10X Genomics clustering and "reclustering" is for reclustering the cells according to the given subset.
#' @param GEX.out.directory.list Character vector that give the path to filtered_feature_bc_matrix data.
#' @param vgm_VDJ Data frame containing cell of interest and x and y coordinates and GEX_barcode.
#' @param vgm_cluster Data frame containing GEX barcode and cluster given by 10X Genomics. Only needed if cluster parameter is set to "GEX_cluster".
#' @param sample_names Character vector containing the name of the sample.
#' @param bcs_merge Data frame containing imagerow, imagecol and barcode of the cells belonging to the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 10.
#' @param images_tibble Tbl-df containing the sample name, grob, height and width of the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 5.
#' @param title Character vector to name the plot.
#' @param size Number, to define the size of the text, default = 15.
#' @param legend_title Character vector to name the legend scale.
#' @return If plotting = TRUE, returns a list containing [[1]] the plot of the selected cells according to their group, [[2]] a data frame that contains the column seurat_clusters with the new cluster. If plotting = FALSE, it returns just the data frame.
#' @export
#' @examples
#' \dontrun{
#' #Clustering of whole cells regardless of cell type
#' GEX_cluster_B_cells<-Spatial_cluster(cluster = "GEX_cluster",
#' vgm_cluster = vgm_with_simulated_VDJ$spatial$cluster[[1]],
#' vgm_VDJ = vgm_with_simulated_VDJ$VDJ,
#' GEX.out.directory.list = GEX.out.directory.list[[1]],images_tibble=scaling_parameters[[5]],
#' bcs_merge=scaling_parameters[[10]], title = "B cells",
#' sample_names = sample_names, legend_title = "GEX clusters" )
#' GEX_cluster_B_cells[[1]]
#'
#' #Reclustering with only B cells
#' reclustering_B_cells<-Spatial_cluster(cluster = "reclustering",
#' vgm_VDJ = vgm_with_simulated_VDJ$VDJ,
#' GEX.out.directory.list = GEX.out.directory.list[[1]],
#' images_tibble=scaling_parameters[[5]],bcs_merge=scaling_parameters[[10]],
#' title = "B cells", sample_names = sample_names, legend_title = "Reclustering")
#' reclustering_B_cells[[1]]
#'}

Spatial_cluster<-function(cluster=c("GEX_cluster","reclustering"),GEX.out.directory.list,vgm_VDJ,vgm_cluster,sample_names,bcs_merge,images_tibble,title,size,legend_title){

  if(missing(cluster)) stop("Please choose between GEX_cluster or reclustering as input method for this function")
  if(missing(sample_names))stop("Please provide sample_names input for this function")
  if(missing(title)){
    title = ""
  }
  if(missing(size)){
    size = 15
  }
  if(missing(legend_title)){
    legend_title = ""
  }
  if(missing(vgm_cluster)) {
    vgm_cluster = NULL
  }
  if(missing(vgm_VDJ)) stop("Please provide vgm_VDJ input for this function")
  if(missing(GEX.out.directory.list)) stop("Please provide GEX.out.directory.list input for this function")

  platypus.version <- "v3"

  x = NULL
  y = NULL
  Cluster = NULL
  grob = NULL
  width = NULL
  height = NULL

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

  if(cluster == "GEX_cluster"){
    vgm_cluster$Barcode<-gsub("-1","",as.character(vgm_cluster$Barcode))
    vgm_cluster$Barcode<-gsub("s1_","",as.character(vgm_cluster$Barcode))
    names(vgm_cluster)[1]<-"barcode" #GEX_barcode
    vgm_VDJ$barcode<-gsub("-1","",as.character(vgm_VDJ$barcode))
    vgm_VDJ$barcode<-gsub("s1_","",as.character(vgm_VDJ$barcode))
    vgm_VDJ_seurat_cluster<-merge(vgm_VDJ,vgm_cluster, by = "barcode")
  } else if (cluster == "reclustering"){
    vgm_cluster = NULL
    #Create a Seurat object
    seurat_object<-Seurat::CreateSeuratObject(Seurat::Read10X(GEX.out.directory.list[[1]]))
    #Subset of Seurat object according to the GEX barcode of the cells of interest
    vgm_VDJ_barcode<-vgm_VDJ$barcode
    vgm_VDJ_barcode<-gsub("-1","",as.character( vgm_VDJ_barcode))
    vgm_VDJ_barcode<-paste0(vgm_VDJ_barcode,"-1")
    #New clustering
    subset_seurat_object<-subset(seurat_object,cells=vgm_VDJ_barcode)
    subset_seurat_object <- Seurat::NormalizeData(subset_seurat_object)
    subset_seurat_object<- Seurat::FindVariableFeatures(subset_seurat_object, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(subset_seurat_object)
    subset_seurat_object <- Seurat::ScaleData(subset_seurat_object, features = all.genes)
    subset_seurat_object <- Seurat::RunPCA(subset_seurat_object, npcs = 50, features = Seurat::VariableFeatures(object = subset_seurat_object))
    subset_seurat_object <- Seurat::FindNeighbors(subset_seurat_object, dims = 1:10)
    subset_seurat_object <- Seurat::FindClusters(subset_seurat_object, resolution = 1)
    #Add new clustering to VDJ
    subset_seurat_object <- as.data.frame(subset_seurat_object@active.ident)
    subset_seurat_object$barcode <- row.names(subset_seurat_object)
    #subset_seurat_object$barcode<-gsub("-1","",as.character(subset_seurat_object$barcode))
    names(subset_seurat_object)[1]<-"Cluster"
    names(subset_seurat_object)[2]<-"barcode"
    subset_seurat_object$barcode<-gsub("-1","",as.character( subset_seurat_object$barcode))
    vgm_VDJ$barcode<-gsub("-1","",as.character(vgm_VDJ$barcode))##
    vgm_VDJ_seurat_cluster<-merge(vgm_VDJ,subset_seurat_object, by = "barcode")
  }

  #Plotting
  p<-ggplot2::ggplot(data = vgm_VDJ_seurat_cluster, ggplot2::aes(x=x,y=y, fill = as.factor(Cluster)))+
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

  return(list(p,vgm_VDJ_seurat_cluster))
}
