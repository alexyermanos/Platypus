#' Plotting the phylogenetic network of a clonotype based on the somatic hypermutations of the immune repertoire sequences on a spatial image.
#' @param simulation Logical operator, to describe which type of data we want to plot, TRUE if the data are output of Echidna simulation and FALSE if the we use real dataset.
#' @param AbForest_output Igraph of phylogenetic tree of a clonotype of interest found in the large list output from AntibodyForest function, only needed if we use real dataset.
#' @param VDJ Data frame containing VDJ information, found in the vgm made by platypus. It must have x and y coordinates column and the column containing the factor to plot.
#' @param nb_clonotype Numeric, value which designates the clonotype we want to study if we use simulated data (Echidna output).
#' @param simulated_VDJ Large list, output of Echidna simulate_repertoire function. Only needed if we use simulated data.
#' @param tracking_type Integer, to define how daughter cells are linked to mother cells.If "all" parameter it means that each daughter cell is link by all these potential mother cells and if "closest" parameter, only closest potential mother cell is link to the daughter cell.
#' @param sample_names Character vector containing the name of the sample.
#' @param bcs_merge Data frame containing imagerow, imagecol and barcode of the cells belonging to the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 10.
#' @param images_tibble Tbl-df containing the sample name, grob, height and width of the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 5.
#' @param title Character vector to name the plot.
#' @param size Number, to define the size of the text, default = 15.
#' @param legend_title Character vector to name the legend scale.
#' @return Plot of phylogenetic network of a clonotype of interest.
#' @export
#' @examples
#' \dontrun{
#'Spatial_evolution_of_clonotype_plot(simulation = FALSE,
#'tracking_type = "closest",AbForest_output = forest$s1$clonotype10,VDJ=vgm$VDJ,
#'sample_names = sample_names, images_tibble = scaling_parameters[[5]],
#'bcs_merge = scaling_parameters[[10]],
#'title = "Tracking evolution of clonotype 10", legend_title = "nb of SHM" )
#'
#'Spatial_evolution_of_clonotype_plot(simulation = FALSE,tracking_type = "all",
#'AbForest_output = forest$s1$clonotype10,VDJ=vgm$VDJ,
#'sample_names = sample_names, images_tibble = scaling_parameters[[5]],
#'bcs_merge = scaling_parameters[[10]],
#'title = "Tracking evolution of clonotype 10", legend_title = "nb of SHM" )
#'
#'Spatial_evolution_of_clonotype_plot(simulation = TRUE,tracking_type = "closest",
#'nb_clonotype = 11 ,simulated_VDJ = simulated_B_cells_VDJ,
#'VDJ =vgm_with_simulated_VDJ$VDJ,bcs_merge = bcs_merge,
#'images_tibble = scaling_parameters[[5]],title = "B cell density",
#'legend_title = "nb_of_SHM",sample_names=sample_names)
#'
#'Spatial_evolution_of_clonotype_plot(simulation = TRUE,tracking_type = "all",
#'nb_clonotype = 11 ,simulated_VDJ = simulated_B_cells_VDJ,
#'VDJ =vgm_with_simulated_VDJ$VDJ,bcs_merge = bcs_merge,
#'images_tibble = scaling_parameters[[5]],title = "B cell density",
#'legend_title = "nb_of_SHM",sample_names=sample_names)
#'}

Spatial_evolution_of_clonotype_plot<-function(simulation =c(TRUE,FALSE),
                                              AbForest_output,
                                              VDJ,nb_clonotype,
                                              simulated_VDJ,
                                              tracking_type = c("closest","all"),
                                              sample_names,
                                              bcs_merge,
                                              images_tibble,
                                              title,
                                              size,
                                              legend_title){

  if(missing(sample_names))stop("Please provide sample_names input for this function")
  if(missing(images_tibble))stop("Please provide images_tibble input for this function")
  if(missing(bcs_merge))stop("Please provide bcs_merge input for this function")
  if(missing(tracking_type)){
    tracking_type = "all"
  }
  if(missing(legend_title)){
    legend_title=""
  }
  if(missing(title)){
    title=""
  }
  if(missing(size)){
    size=15
  }

  platypus.version <- "v3"

  from = NULL
  to = NULL
  x = NULL
  y = NULL
  grob = NULL
  width = NULL
  height = NULL
  Daughter = NULL
  mother_label = NULL
  x0 = NULL
  y0 = NULL
  x1 = NULL
  y1 = NULL

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

  #Real data--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(simulation == FALSE){
    #Set up parameter only needed when simulation = FALSE and put to NULL parameter only needed if simulation = TRUE
    if(missing(AbForest_output))stop("Please provide AbForest_output input for this function")
    if(missing(VDJ))stop("Please provide VDJ input for this function")
    nb_clonotype=NULL
    simulated_VDJ=NULL
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
    available_cells<-merge(VDJ, tree_dataframe, by = "barcode")
    #Delete the from germline label
    germline_label<-dplyr::filter(vertices_dataframe, vertices_dataframe$node_type == "germline")
    germline_label<-germline_label$label
    tree_pathway<-dplyr::filter(tree_pathway, from != germline_label)
    tree_pathway<-dplyr::filter(tree_pathway, to != germline_label)
    if(length(tree_pathway$from)==0){
      plot<-ggplot2::ggplot(data =available_cells, ggplot2::aes(x=x,y=y, fill = as.factor(label)))+
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
      stop(return(plot))
    }
    #Mother cell is the closest potential mother cell of the daughter cell
    if(tracking_type == "closest"){
      #Loop to calculate euclidan matrix and link mother and daughter cells
      final_arrow_dataframe<-list()
      for (l in 1:length(unique(tree_pathway$from))) {
        subset_available_cells<-list()
        from_cell<-tree_pathway$from[[1]]

        #remove form tree_pathway
        tree_law<-subset(tree_pathway, from ==from_cell )
        rownumber<-as.numeric(rownames(tree_law))
        tree_pathway<-tree_pathway[-c(rownumber),]
        to_cell<-tree_law$to

        #Assigne mother or daughter for each selected cells
        subset_available_cells <- available_cells[available_cells$label %in% to_cell | available_cells$label%in%from_cell,]
        nb_daughter_cell<-0
        for (i in 1:length(subset_available_cells$barcode)){
          if(subset_available_cells$label[[i]] == from_cell){
            subset_available_cells$mother_or_daughter[[i]]<-"mother"
          }else{
            subset_available_cells$mother_or_daughter[[i]]<-"daughter"
            nb_daughter_cell<-nb_daughter_cell+1
          }
        }
        #Link a daughter cell to a mother cell
        mother_cell_final<-list()
        daughter_cell_final<-list()
        for(k in 1:nb_daughter_cell){
          #Extract coordinates for euclidan matrix
          mother_cells<-dplyr::filter(subset_available_cells,subset_available_cells$mother_or_daughter == "mother")
          mother_cells<-dplyr::select(mother_cells, barcode,x,y)
          daughter_cells<-dplyr::filter(subset_available_cells,subset_available_cells$mother_or_daughter == "daughter")
          daughter_cells<-dplyr::select(daughter_cells, barcode, x,y)

          #Matrix euclidean to find the closest cell
          euclidian_matrix<-matrix(data=NA, nrow = length(mother_cells$barcode), ncol = length(daughter_cells$barcode))
          CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))
          for(i in 1:length(mother_cells$barcode)){
            for(j in 1:length(daughter_cells$barcode)){
              m_x<-mother_cells$x[[i]]
              m_y<-mother_cells$y[[i]]
              m_v<-c(m_x,m_y)
              d_x<-daughter_cells$x[[j]]
              d_y<-daughter_cells$y[[j]]
              d_v<-c(d_x,d_y)
              euclidian_matrix[i,j]<-CalculateEuclideanDistance(m_v,d_v)
            }
          }
          rownames(euclidian_matrix)<-mother_cells$barcode
          colnames(euclidian_matrix)<-daughter_cells$barcode

          #Extract mother and daughter cell from nearest cells
          min_value<-as.data.frame(which(euclidian_matrix == min(euclidian_matrix), arr.ind=TRUE))
          mother_cell_final[[k]]<-rownames(euclidian_matrix)[min_value$row]
          daughter_cell_final[[k]]<-colnames(euclidian_matrix)[min_value$col]

          #Remove daughter cell from subset available cells
          subset_available_cells<- subset(subset_available_cells, barcode != daughter_cell_final[[k]])
        }
        #Results of linking mother and daughter cells
        mother <- t(as.data.frame(mother_cell_final))
        daughter <- t(as.data.frame(daughter_cell_final))
        final <- data.frame(mother, daughter)
        final_arrow_dataframe<-rbind(final_arrow_dataframe,final)
      }
      #Add coordinates for each cells
      names(final_arrow_dataframe)<-c("barcode","Daughter")
      arrow_coordinates<-merge(final_arrow_dataframe, available_cells, by = "barcode")
      arrow_coordinates<-dplyr::select(arrow_coordinates, barcode, x,y,label, Daughter)
      names(arrow_coordinates)<-c("mother","x0","y0","mother_label","barcode")
      arrow_coordinates<-merge(arrow_coordinates,available_cells)
      arrow_coordinates<-dplyr::select(arrow_coordinates,mother,x0,y0,mother_label,barcode,x,y,label)
      names(arrow_coordinates)<-c("mother_barcode","x0","y0","mother_label","daughter_barcode","x1","y1","daughter_label")
    }
    #All mother cells link all daughter cells
    if(tracking_type == "all"){
      #Loop for assignment mother-daughter cells
      final_arrow_dataframe<-list()
      mother_barcode<-list()
      daughter_barcode<-list()
      for (l in 1:length(unique(tree_pathway$from))) {
        subset_available_cells<-list()
        from_cell<-tree_pathway$from[[1]]
        #remove form tree_pathway
        tree_law<-subset(tree_pathway, from ==from_cell )
        to_cell<-tree_law$to
        tree_pathway<-dplyr::filter(tree_pathway, from != from_cell)
        #Assign mother or daughter for each selected cells
        subset_available_cells <- available_cells[available_cells$label %in% to_cell | available_cells$label%in%from_cell,]
        nb_daughter_cell<-0
        for (i in 1:length(subset_available_cells$barcode)){
          if(subset_available_cells$label[[i]] == from_cell){
            subset_available_cells$mother_or_daughter[[i]]<-"mother"
          }else{
            subset_available_cells$mother_or_daughter[[i]]<-"daughter"
            nb_daughter_cell<-nb_daughter_cell+1
          }
        }
        #Link a daughter cell to a mother cell
        mother_cells<-dplyr::filter(subset_available_cells,subset_available_cells$mother_or_daughter == "mother")
        mother_cells<-dplyr::select(mother_cells, barcode,x,y)
        daughter_cells<-dplyr::filter(subset_available_cells,subset_available_cells$mother_or_daughter == "daughter")
        daughter_cells<-dplyr::select(daughter_cells, barcode, x,y)

        for(k in 1:length(mother_cells$barcode)){
          for(m in 1:length(daughter_cells$barcode)){
            mother_barcode[[m]]<-mother_cells$barcode[[k]]
            daughter_barcode[[m]]<-daughter_cells$barcode[[m]]
          }
          mother <- t(as.data.frame(mother_barcode))
          daughter <- t(as.data.frame(daughter_barcode))
          fdata<-data.frame(mother,daughter)
          final_arrow_dataframe<-rbind(final_arrow_dataframe,fdata)
          mother_barcode<-list()
          daughter_barcode<-list()
        }
      }

      #Add coordinates for each cells
      names(final_arrow_dataframe)<-c("barcode","Daughter")
      arrow_coordinates<-merge(final_arrow_dataframe, available_cells, by = "barcode")
      arrow_coordinates<-dplyr::select(arrow_coordinates, barcode, x,y,label, Daughter)
      names(arrow_coordinates)<-c("mother","x0","y0","mother_label","barcode")
      arrow_coordinates<-merge(arrow_coordinates,available_cells)
      arrow_coordinates<-dplyr::select(arrow_coordinates,mother,x0,y0,mother_label,barcode,x,y,label)
      names(arrow_coordinates)<-c("mother_barcode","x0","y0","mother_label","daughter_barcode","x1","y1","daughter_label")
    }

    #Plotting
    plot<-ggplot2::ggplot(mapping = ggplot2::aes(x,y))+
      geom_spatial(data=images_tibble[1,], ggplot2::aes(grob=grob), x=0.5, y=0.5)+ #add spatial image
      ggplot2::xlim(0,max(bcs_merge %>%
                            dplyr::filter(sample ==sample_names[1]) %>%
                            dplyr::select(width)))+
      ggplot2::ylim(max(bcs_merge %>%
                          dplyr::filter(sample ==sample_names[1]) %>%
                          dplyr::select(height)),0)+
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::ggtitle(sample_names[1], title)+
      ggplot2::theme(axis.text=ggplot2::element_text(size=15),
            axis.title=ggplot2::element_text(size=15))+
      ggplot2::geom_point(data = available_cells, ggplot2::aes(x,y,fill=factor(distance_from_germline)),shape = 21, colour = "black", size = 1.5, stroke = 0.5)+ #cells with 0 SHM without germline
      ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(reverse=TRUE))+
      ggplot2::geom_segment(data = arrow_coordinates, ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1), arrow = ggplot2::arrow(length = grid::unit(0.01, "npc")),colour = "#FFFFFF")+
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size=3)))+
      ggplot2::theme_set(ggplot2::theme_bw(base_size = 10))+
      ggplot2::labs(fill = "Number of SHM")+
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank())
    return(plot)
  }
  #Simulated data with Echidna--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(simulation == TRUE){
    #Setup parameter
    if(missing(nb_clonotype))stop("Please provide nb_clonotype input for this function")
    if(missing(simulated_VDJ))stop("Please provide simulated_VDJ input for this function")

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
    VDJ$orig_barcode<-paste0(VDJ$orig_barcode,"-1")
    available_cells<-merge(VDJ,history, by = "orig_barcode")
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
    #3)Tree with only nodes with germline node relationship but germline node empty
    if(length(tree_pathway$from)==0){
      plot<-ggplot2::ggplot(data =available_cells, ggplot2::aes(x=x,y=y, fill = as.factor(label)))+
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
      stop(return(plot))
    }
    #4)all parameter
    if(tracking_type == "all"){
      final_arrow_dataframe<-list()
      mother_barcode<-list()
      daughter_barcode<-list()
      for (l in 1:length(unique(tree_pathway$from))) {
        subset_available_cells<-list()
        from_cell<-tree_pathway$from[[1]]
        #remove form tree_pathway
        tree_law<-subset(tree_pathway, from ==from_cell )
        tree_pathway<-dplyr::filter(tree_pathway,from!=from_cell)
        to_cell<-tree_law$to
        #Assigne mother or daughter for each selected cells
        subset_available_cells <- available_cells[available_cells$label %in% to_cell | available_cells$label%in%from_cell,]
        nb_daughter_cell<-0
        for (i in 1:length(subset_available_cells$barcode)){
          if(subset_available_cells$label[[i]] == from_cell){
            subset_available_cells$mother_or_daughter[[i]]<-"mother"
          }else{
            subset_available_cells$mother_or_daughter[[i]]<-"daughter"
            nb_daughter_cell<-nb_daughter_cell+1
          }
        }
        #Link a daughter cell to a mother cell
        mother_cells<-dplyr::filter(subset_available_cells,subset_available_cells$mother_or_daughter == "mother")
        mother_cells<-dplyr::select(mother_cells, barcode,x,y)
        daughter_cells<-dplyr::filter(subset_available_cells,subset_available_cells$mother_or_daughter == "daughter")
        daughter_cells<-dplyr::select(daughter_cells, barcode, x,y)

        for(k in 1:length(mother_cells$barcode)){
          for(m in 1:length(daughter_cells$barcode)){
            mother_barcode[[m]]<-mother_cells$barcode[[k]]
            daughter_barcode[[m]]<-daughter_cells$barcode[[m]]
          }
          mother <- t(as.data.frame(mother_barcode))
          daughter <- t(as.data.frame(daughter_barcode))
          fdata<-data.frame(mother,daughter)
          final_arrow_dataframe<-rbind(final_arrow_dataframe,fdata)
          mother_barcode<-list()
          daughter_barcode<-list()
        }
      }

      #3)Add coordinates for each cells
      names(final_arrow_dataframe)<-c("barcode","Daughter")
      arrow_coordinates<-merge(final_arrow_dataframe, available_cells, by = "barcode")
      arrow_coordinates<-dplyr::select(arrow_coordinates, barcode, x,y,label, Daughter)
      names(arrow_coordinates)<-c("mother","x0","y0","mother_label","barcode")
      arrow_coordinates<-merge(arrow_coordinates,available_cells)
      arrow_coordinates<-dplyr::select(arrow_coordinates,mother,x0,y0,mother_label,barcode,x,y,label)
      names(arrow_coordinates)<-c("mother_barcode","x0","y0","mother_label","daughter_barcode","x1","y1","daughter_label")
    }
    #5) Closest
    if(tracking_type == "closest"){
      final_arrow_dataframe<-list()

      for (l in 1:length(unique(tree_pathway$from))) {
        subset_available_cells<-list()
        from_cell<-tree_pathway$from[[1]]
        #remove form tree_pathway
        tree_law<-subset(tree_pathway, from ==from_cell )
        tree_pathway<-dplyr::filter(tree_pathway, from != from_cell)
        to_cell<-tree_law$to
        #Assigne mother or daughter for each selected cells
        subset_available_cells <- available_cells[available_cells$label %in% to_cell | available_cells$label%in%from_cell,]
        nb_daughter_cell<-0
        for (i in 1:length(subset_available_cells$barcode)){
          if(subset_available_cells$label[[i]] == from_cell){
            subset_available_cells$mother_or_daughter[[i]]<-"mother"
          }else{
            subset_available_cells$mother_or_daughter[[i]]<-"daughter"
            nb_daughter_cell<-nb_daughter_cell+1
          }
        }
        #Link a daughter cell to a mother cell
        #nb_daughter_cell<-length(subset_available_cells$mother_or_daughter == "daughter")
        mother_cell_final<-list()
        daughter_cell_final<-list()
        for(k in 1:nb_daughter_cell){
          #Extract coordinates for euclidan matrix
          mother_cells<-dplyr::filter(subset_available_cells,subset_available_cells$mother_or_daughter == "mother")
          mother_cells<-dplyr::select(mother_cells, barcode,x,y)
          daughter_cells<-dplyr::filter(subset_available_cells,subset_available_cells$mother_or_daughter == "daughter")
          daughter_cells<-dplyr::select(daughter_cells, barcode, x,y)

          #Matrix euclidean to find the closest cell
          euclidian_matrix<-matrix(data=NA, nrow = length(mother_cells$barcode), ncol = length(daughter_cells$barcode))
          CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))
          for(i in 1:length(mother_cells$barcode)){
            for(j in 1:length(daughter_cells$barcode)){
              m_x<-mother_cells$x[[i]]
              m_y<-mother_cells$y[[i]]
              m_v<-c(m_x,m_y)
              d_x<-daughter_cells$x[[j]]
              d_y<-daughter_cells$y[[j]]
              d_v<-c(d_x,d_y)
              euclidian_matrix[i,j]<-CalculateEuclideanDistance(m_v,d_v)
            }
          }
          rownames(euclidian_matrix)<-mother_cells$barcode
          colnames(euclidian_matrix)<-daughter_cells$barcode

          #Extract mother and daughter cell from nearest cells
          min_value<-as.data.frame(which(euclidian_matrix == min(euclidian_matrix), arr.ind=TRUE))
          mother_cell_final[[k]]<-rownames(euclidian_matrix)[min_value$row]
          daughter_cell_final[[k]]<-colnames(euclidian_matrix)[min_value$col]

          #Remove daughter cell from subset available cells
          subset_available_cells<- subset(subset_available_cells, barcode != daughter_cell_final[[k]])
        }
        #Results of linking mother and daughter cells
        mother <- t(as.data.frame(mother_cell_final))
        daughter <- t(as.data.frame(daughter_cell_final))
        final <- data.frame(mother, daughter)
        final_arrow_dataframe<-rbind(final_arrow_dataframe,final)
      }

      #3)Add coordinates for each cells
      names(final_arrow_dataframe)<-c("barcode","Daughter")
      arrow_coordinates<-merge(final_arrow_dataframe, available_cells, by = "barcode")
      arrow_coordinates<-dplyr::select(arrow_coordinates, barcode, x,y,label, Daughter)
      names(arrow_coordinates)<-c("mother","x0","y0","mother_label","barcode")
      arrow_coordinates<-merge(arrow_coordinates,available_cells)
      arrow_coordinates<-dplyr::select(arrow_coordinates,mother,x0,y0,mother_label,barcode,x,y,label)
      names(arrow_coordinates)<-c("mother_barcode","x0","y0","mother_label","daughter_barcode","x1","y1","daughter_label")
    }

    #Plotting
    plot<-ggplot2::ggplot(mapping = ggplot2::aes(x,y))+
      geom_spatial(data=images_tibble[1,], ggplot2::aes(grob=grob), x=0.5, y=0.5)+ #add spatial image
      ggplot2::xlim(0,max(bcs_merge %>%
                            dplyr::filter(sample ==sample_names[1]) %>%
                            dplyr::select(width)))+
      ggplot2::ylim(max(bcs_merge %>%
                          dplyr::filter(sample ==sample_names[1]) %>%
                          dplyr::select(height)),0)+
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::ggtitle(sample_names[1], title)+
      ggplot2::theme(axis.text=ggplot2::element_text(size=15),
            axis.title=ggplot2::element_text(size=15))+
      ggplot2::geom_point(data = available_cells, ggplot2::aes(x,y,fill=factor(distance_from_germline)),shape = 21, colour = "black", size = 3, stroke = 0.5)+ #cells with 0 SHM without germline
      ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(reverse=TRUE))+
      ggplot2::geom_segment(data = arrow_coordinates, ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1), arrow = ggplot2::arrow(length = grid::unit(0.01, "npc")),colour = "#FFFFFF")+
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size=3)))+
      ggplot2::theme_set(ggplot2::theme_bw(base_size = 10))+
      ggplot2::labs(fill = legend_title)+
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank())
    return(plot)
  }
}
