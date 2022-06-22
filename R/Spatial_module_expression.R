#' Plotting the expression of a gene module.
#' @param sample_names Character vector containing the name of the sample.
#' @param gene.set Charcter vector containing the markers name.
#' @param GEX.out.directory.list Character vector that give the path to filtered_feature_bc_matrix data.
#' @param bcs_merge Data frame containing imagerow, imagecol and barcode of the cells belonging to the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 10.
#' @param images_tibble Tbl-df containing the sample name, grob, height and width of the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 5.
#' @param title Character vector to name the plot.
#' @param size Number, to define the size of the text, default = 15.
#' @param threshold Number, to define the threshold. If threshold = No, plot of the module and if threshold is a number, plot show the cells above the threshold.
#' @param legend_title Character vector to name the legend scale.
#' @return Returns a plot gene module expression.
#' @export
#' @examples
#' \dontrun{
#' gene.set <- list() # make empty list 
#' gene.set[[1]] <- c("CD19","XBP1","SDC1") # put gene set in list
#' 
#' #Without expression threshold
#' Spatial_module_expression(sample_names = sample_names,gene.set = gene.set, 
#' GEX.out.directory.list = GEX.out.directory.list[[1]],bcs_merge = scaling_parameters[[10]],
#' images_tibble = scaling_parameters[[5]], threshold = "No")
#' 
#' #With expression threshold
#' Spatial_module_expression(sample_names = sample_names,gene.set = gene.set, 
#' GEX.out.directory.list = GEX.out.directory.list[[1]],bcs_merge = scaling_parameters[[10]],
#' images_tibble = scaling_parameters[[5]], threshold = 1)
#' }

Spatial_module_expression<-function(sample_names,gene.set,GEX.out.directory.list,bcs_merge,images_tibble,title,size,threshold,legend_title){
  
  seurat_object<-NULL
  seurat_data_frame<-NULL
  coordinates<-NULL
  module_data_frame<-NULL
  
  if (missing(threshold)){
    threshold = "No"
  }
  if (missing(title)){
    title = ""
  }
  if (missing(legend_title)){
    legend_title = ""
  }
  if(missing(size)){
    size = 15
  }
  if(missing(sample_names)) stop("Please provide sample_names input for this function")
  if(missing(images_tibble)) stop("Please provide images_tibble input for this function")
  if(missing(bcs_merge)) stop("Please provide bcs_merge input for this function")
  if(missing(gene.set)) stop("Please provide gene.set input for this function")
  if(missing(GEX.out.directory.list)) stop("Please provide GEX.out.directory.list input for this function")
  
  platypus.version <- "v3" 
  
  #Seurat object
  seurat_object<-CreateSeuratObject(Read10X(GEX.out.directory.list))
  #this Adds the module score to your object as a Feature 
  seurat_object <- AddModuleScore(object = seurat_object,features = list(gene.set[[1]]))
  seurat_data_frame<-seurat_object@meta.data
  seurat_data_frame$barcode <- row.names(seurat_data_frame)
  names(seurat_data_frame)[4]<-"module"
  seurat_data_frame<-select(seurat_data_frame, barcode, module)
  coordinates<-bcs_merge
  coordinates<-select(coordinates, imagerow, imagecol, barcode)
  module_data_frame<-merge(seurat_data_frame, coordinates, by="barcode")
  
  if (threshold == "No"){
    #Colors
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    #Plotting
    p<-ggplot(data = module_data_frame, aes(x=imagecol,y=imagerow, fill = module))+
      geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape=21, colour = "black", size = 1.75, stroke = 0.5)+
      coord_cartesian(expand=FALSE)+
      scale_fill_gradientn(colours = myPalette(100))+
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
      labs(fill = legend_title)+
      theme_set(theme_bw(base_size = size))+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank())
  } else {
    for (i in 1:length(module_data_frame$barcode)){
      if (module_data_frame$module[[i]]>=threshold){
        module_data_frame$threshold[[i]] = 1
      } else{
        module_data_frame$threshold[[i]] = 0
      }
    }
    module_data_frame$module<-as.numeric(module_data_frame$threshold)
    module_data_frame<-filter(module_data_frame,module==1)
    p<-ggplot(data = module_data_frame, aes(x=imagecol,y=imagerow, fill = as.factor(module)))+
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
  return(p)
}
