#' Plotting number of somatic hypermutation of clones compare to the germline sequence of the clonotype.
#' @param simulation Logical operator, to describe which type of data we want to plot, TRUE if the data are output of Echidna simulation and FALSE if the we use real dataset.
#' @param vgm_VDJ Data frame containing cell of interest and x and y coordinates and GEX_barcode.
#' @param AbForest_output Igraph of phylogenetic tree of a clonotype of interest found in the large list output from AntibodyForest function, only needed if we use real dataset.
#' @param nb_clonotype Numeric, value which designates the clonotype we want to study if we use simulated data (Echidna output).
#' @param simulated_VDJ Large list, output of Echidna simulate_repertoire function. Only needed if we use simulated data.
#' @param sample_names Character vector containing the name of the sample.
#' @param bcs_merge Data frame containing imagerow, imagecol and barcode of the cells belonging to the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 10.
#' @param images_tibble Tbl-df containing the sample name, grob, height and width of the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 5.
#' @param title Character vector to name the plot.
#' @param size Number, to define the size of the text, default = 15.
#' @param legend_title Character vector to name the legend scale.
#' @return Spatial plot with cells colored by number of somatic hypermutation
#' @export
#' @examples
#' \dontrun{
#' Spatial_nb_SHM_compare_to_germline_plot(simulation = FALSE,AbForest_output=forest[[1]][[2]], vgm_VDJ = vgm$VDJ,
#' images_tibble = scaling_parameters[[5]],bcs_merge = scaling_parameters[[10]],sample_names = sample_names,
#' title = "Number of SHM of clonotype 10", legend_title = "nb of SHM")
#'}

Spatial_nb_SHM_compare_to_germline_plot<-function(simulation = c(TRUE,FALSE),vgm_VDJ,AbForest_output,nb_clonotype,simulated_VDJ,sample_names,bcs_merge,images_tibble,title,size,legend_title){
  if(missing(vgm_VDJ)) stop("Please provide vgm_VDJ input for this function")
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
  if(missing(simulation)){
    simulation = FALSE
  }
  
  platypus.version <- "v3" 
  
  from = NULL
  x = NULL
  y = NULL
  grop = NULL
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
  
  #Simulated Data with Echidna------------------------------------------------------------------------------------------------------------------
  if(simulation == TRUE){
    #Set to null parameter needed only if simulation == FALSE
    AbForest_output=NULL
    #Set up parameter only needed if simulation == TRUE
    if(missing(simulated_VDJ)) stop("Please provide the simulate VDJ, output of Echidna simulation for this function")
    if(missing(nb_clonotype)) stop("Please provide the clonotype_nb input for this function")
    #Formation of tree path
    clonotype_id<-nb_clonotype
    clonotype_id<-paste0("clonotype",nb_clonotype)
    
    clonoselect<-simulated_VDJ$all_contig_annotations$raw_clonotype_id==clonotype_id
    
    barcode<-simulated_VDJ$all_contig_annotations$barcode[clonoselect]
    sequence<-simulated_VDJ$all_contig$seq[clonoselect]
    
    history<-simulated_VDJ$history[simulated_VDJ$history$barcode.history %in% barcode,]
    
    colnames(history)<-c("label","orig_barcode","network_sequences","distance_from_germline","trans_state_his" )
    .RM.EMPTY.NODE<-function(No,size0){
      if(length(size0)>0){
        size0<-rev(size0)
        empty_del<-as.data.frame(t(matrix(No,2)))
        for (e in size0){
          mom<-empty_del[,1][which(empty_del[,2]==e)]
          empty_del<-empty_del[-which(empty_del[,2]==e),]
          empty_del[,1][which(empty_del[,1]==e)]<-mom
          empty_del[empty_del>e]<-empty_del[empty_del>e]-1
        }
        No<-as.vector(t(empty_del))
      }
      
      return(No)
    }
    
    igraph_index_rm_empty<-function(igraph.index,history){
      #igraph
      size.of.vertex<-as.data.frame(stats::aggregate(barcode.history~seq.number,history,length,na.action = stats::na.pass))
      #NA count as 1
      size.of.vertex1<-as.data.frame(stats::aggregate(barcode.history~seq.number,history,length))
      size.of.vertex<-merge(size.of.vertex1,size.of.vertex,all=T, by="seq.number")[,-3]
      size.of.vertex[!(size.of.vertex$seq.number %in% size.of.vertex1$seq.number),2]<-0
      #NA count as 0
      colnames(size.of.vertex)<-c("seq.number","size")
      
      igraph.index.attr<-data.frame(igraph.index)
      colnames(igraph.index.attr)<-"seq.number"
      igraph.index.attr<-subset(igraph.index.attr,!duplicated(igraph.index.attr$seq.number))
      #size of vertex
      igraph.index.attr<-merge(size.of.vertex,igraph.index.attr,by="seq.number",all.y=T)
      igraph.index.attr$size[1]<-0
      
      
      igraph.index.jr<-data.frame(igraph.index)
      colnames(igraph.index.jr)<-"seq.number"
      igraph.index.jr$No<-match(x=igraph.index.jr$seq.number,table=igraph.index.attr$seq.number)
      
      #find empty
      size0<-which(igraph.index.attr$size==0)
      
      No<-.RM.EMPTY.NODE(igraph.index.jr$No,size0[-1])
      
      table<-cbind(igraph.index.attr[igraph.index.attr$size!=0,],No=2:(sum(igraph.index.attr$size!=0)+1))
      table[nrow(table)+1,]<-c(0,0,1)
      output<-table$seq.number[ match(No,table$No)]
      
      
      return(output)
    }
    igraph_index<-igraph_index_rm_empty(igraph.index = simulated_VDJ$igraph.index[[nb_clonotype]],history = simulated_VDJ$history)
    tree_pathway<-as.data.frame(t(matrix(igraph_index,2)))
    colnames(tree_pathway)<-c("from","to")
    barcode<-unique(barcode)
    vgm_VDJ$orig_barcode<-paste0(vgm_VDJ$orig_barcode,"-1")
    available_cells<-merge(vgm_VDJ,history, by = "orig_barcode")
    history$orig_barcode<-available_cells$barcode
    names(available_cells)[2]<-"old_barcode"
    names(available_cells)[1]<-"barcode"
    #2)modify Tree_pathway
    germline<-any(available_cells$label ==0)
    if(germline == FALSE ){
      germline_nb<-dplyr::filter(tree_pathway, from == 0)
      germline_nb<-germline_nb$to
      germline_seq<-dplyr::filter(available_cells, label==germline_nb)
      germline_seq<-unique(germline_seq$network_sequences)
      tree_pathway<-dplyr::filter(tree_pathway, tree_pathway$from !=0)
    }
    if(germline == TRUE){
      tree_pathway<-tree_pathway
      germline_seq<-dplyr::filter(available_cells, label == 0)
      germline_seq<-unique(germline_seq$VDJ_sequence_nt_raw)
    }
    #Add distance_from_germline
    for(i in 1:length(available_cells$barcode)){
      available_cells$distance_from_germline[[i]]<-get.seq.distance(germline=germline_seq, sequence=available_cells$network_sequences[[i]])
    }
    #Plotting
    plot<-ggplot2::ggplot(data = available_cells, ggplot2::aes(x=x,y=y, fill = as.factor(distance_from_germline)))+
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
  }
  #Real data----------------------------------------------------------------------------------------------------------------------------------------------
  if(simulation == FALSE){
    #Set to null parameter needed only if simulation == TRUE
    simulated_VDJ = NULL
    #Set up parameter only needed if simulation == FALSE
    if(missing(AbForest_output)) stop("Please provide AbForest_output input for this function")
    #Preparation of VDJ based on AntibodyForest output containing tree_path and vertices_dataframe 
    vertices_dataframe<-igraph::as_data_frame(AbForest_output, what = "vertices")
    tree_pathway<-igraph::as_data_frame(AbForest_output, what = "edges")
    tree_dataframe<-list()
    for (i in 1:length(vertices_dataframe$cell_barcodes)){
      barcode<-as.data.frame(vertices_dataframe$cell_barcodes[[i]])
      names(barcode)<-"barcode"
      l<-vertices_dataframe$label[[i]]
      label<-as.data.frame(rep(l,length(vertices_dataframe$cell_barcodes[[i]])))
      names(label)<-"label"
      m<-vertices_dataframe$distance_from_germline[[i]]
      distance_from_germline<-as.data.frame(rep(m,length(vertices_dataframe$cell_barcodes[[i]])))
      names(distance_from_germline)<-"distance_from_germline"
      dataframe<-cbind(barcode, label,distance_from_germline)
      tree_dataframe<-rbind(tree_dataframe, dataframe)
    }
    new_VDJ<-merge(vgm_VDJ, tree_dataframe, by = "barcode")
    #Plotting
    plot<-ggplot2::ggplot(data = new_VDJ, ggplot2::aes(x=x,y=y, fill = as.factor(distance_from_germline)))+
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
  }
  return(plot)
}







