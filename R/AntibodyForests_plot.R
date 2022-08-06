#' Custom plots for trees/networks created with AntibodyForests


#'@description AntibodyForests_plot takes the input of AntibodyForests and outputs a list of plot-ready graphs (inside the AntibodyForests object) to be further used with plot(). Plots can also be automicatically saved to pdf via the save.pdf parameter. The resulting igraph object have their node/edge colors/shapes/sizes added following the specific parameters in the AntibodyForests_plot function.

#' @param network.list nested list of igraph objects, as obtained from the AntibodyForests function. input[[1]][[2]] represents the igraph object for the first sample, second clonotype, if the normal AntibodyForests parameters are used.
#' @param graph.type string - graph to be plotted - 'tree' for the normal tree plots, 'heterogeneous' for the single-cell grahs (call AntibodyForests_heterogeneous before), 'phylo' to plot the phylo objects (call AntibodyForests_phylo before), 'dynamic' for the dynamic/temporal graphs.
#' @param node.color string specifying the name of the igraph vertex attribute (or the original vgm[[1]]/VDJ column name used in the node.features parameter of AntibodyForests) to be used for coloring the nodes.
#' If the node.shape is also 'pie', the resulting nodes will be a pie chart of the per-node values denoted by the node.color parameter (if there are multiple different values per node - e.g., a node denotes a sequence, which can be further traced back to multiple cells with different barcodes and with potentially different transcriptomic clusters - multiple feature values per node).
#' If the node.color parameter is NULL, the default node color will be '#FFCC00' (yellow) for sequence nodes, lighter gray for intermediate/inferred nodes, and darker gray for germline nodes.
#' @param node.label string - 'cells' to label the nodes by the number of cells with that specific sequence, 'rank' for cell count-based ranking of the nodes. If NULL, will not add number labels to the sequence nodes.
#' @param node.shape string specifying the the name of the igraph vertex attribute (or the original vgm[[1]]/VDJ column name used in the node.features parameter of AntibodyForests) to be used for node shapes. Shapes will be assigned per unique value from the values of this column/vertex attribute. There is a maximum of 7 unique shapes to be chosen from, therefore node.shapes should be used for features with less than 7 unique values.
#' If the node.shape parameter is null, the node shape for all nodes will default to 'circle'
#' @param node.size string denoting either a specific method of scaling the node sizes of the input graph or a specific vertex attribute of numeric values (added via the node.features parameter of the AntibodyForests from the original vgm[[1]]/VDJ dataframe) to be used for node sizes.
#' If NULL, then the sizes will be equal to node.scale.factor * 1. If scaleByEigen, node sizes will be scaled by the eigenvector centrality of each node; scaleByCloseness - closeness centrality; scaleByBetweenness - betweenness centrality; scaleByExpansion - by the sequence frequency of each node, as originally calculated by the AntibodyForests function.
#' @param max.node.size Maximum size of any given node
#' @param node.scale.factor integer to further refine the size of each node. Each vertex size will be multiplied by scaling factor.
#' @param edge.length either NULL - edge lengths are constant, scaleByWeight - edge lengths will be scaled by the edge weight attribute (the string distance between pairs of sequences/nodes, as calculated by the AntibodyForests function) - larger weights/string distances = longer edges/nodes further apart, or a specific igraph edge attribute name.
#' @param edge.width either NULL - edge widths are constant, scaleByWeight - edge widths will be scaled by the edge weight attribute (the string distance between pairs of sequences/nodes, as calculated by the AntibodyForests function) - larger weights/string distances = thinner edges), or a specific igraph edge attribute name.
#' @param path.list named list of igraph paths, as obtained from the AntibodyForests_metrics function.
#' @param specific.node.colors named list of colors to be used for the node.color parameter. If NULL, colors will be automatically added for each unique value.
#' For example, if specific.node.colors=list('yes'='blue', 'no'='red'), then the nodes labeled as 'yes' will be colored blue, the others red.
#' @param specific.node.shapes named list of node shapes to be used for the node.shapes parameter. Must be shapes compatible with igraph objects, use setdiff(igraph::shapes(), "") to get a list of possible values.
#' For example, if specific.node.shapes=list('yes'='circle', 'no'='square'), then the nodes labele as 'yes' will be circles, the others squares.
#' @param specific.edge.colors named list of edge colors to be used for the edge.colors parameter. The names should be the path metrics names obtained from the nested paths list, from AntibodyForests_metrics.
#' For example, if specific.edge.colors=list('longest.path.weighted'='blue', 'shortest.path.unweighted'='red'), the longest weighted paths obtained from AbtibodyForests will be colored blue for each igraph object, the rest will be red.
#' @param color.by.majority boolean - if T, will color the entire network (all nodes) by the dominant/most frequent node feature, as specified in the node.color parameter.
#' @param cell.color string - cell feature column denoting the cell colors - as denoted by the node.features parameter when calling AntibodyForests_heterogeneous.
#' @param specific.cell.colors named list of cell colors and their features (e.g., for Seurat clusters: list(1 = 'red', 2 = 'blue')). Optional (will auto search for unique colors per feature).
#' @param cell.size integer denoting the size of the cell nodes.
#' @param network.layout either NULL - will default to the automatic igraph::layout_nicely(), 'fr' - igraph::layout_with_fr() - for fully connected graphs/graphs with defined connected components, 'tree' - for tree graphs, as obtained from AntibodyForests with network.algorithm='tree'.
#' @param save.pdf boolean - if T, plots will be automatically saved to pdf, in the current working directory; F - normal output of the function (plot-ready igraph object and specific layout). New folders will be created for each sample of the input nested list of igraph objects.
#' @param save.dir path to the directory oin which the network PDFs will be saved.
#' @param show.legend boolean - whether the legend should be showed in the resulting plots
#' @param color.gradient string - defualt: NULL - the feature whose value will be plot as a color gradient
#' @return nested list of plot-ready AntibodyForests objects. Can also save the plots as a PDF file.
#' @export
#' @seealso AntibodyForests, AntibodyForests_metrics
#' @examples
#' \dontrun{
#' AntibodyForests_plot(graphs, node.color='clonotype_id',
#' node.size='scaleByExpansion', network.layout='tree',
#' save.pdf=T)
#'}



AntibodyForests_plot <- function(network.list,
                                 graph.type,
                                 node.color,
                                 node.label,
                                 node.shape,
                                 node.size,
                                 max.node.size,
                                 node.scale.factor,
                                 edge.length,
                                 edge.width,
                                 path.list,
                                 specific.node.colors,
                                 specific.node.shapes,
                                 specific.edge.colors,
                                 color.by.majority,
                                 cell.color,
                                 specific.cell.colors,
                                 cell.size,
                                 network.layout,
                                 save.pdf,
                                 save.dir,
                                 show.legend,
                                 color.gradient
                                 ){

  if(missing(network.list)) stop('Please input a nested list of networks and their corresponding network dataframe, output of AntibodyForests_parallel')
  if(missing(graph.type)) graph.type <- 'tree'
  if(missing(node.color)) node.color <- NULL
  if(missing(node.label)) node.label <- 'cells'
  if(missing(node.shape)) node.shape <- NULL
  if(missing(node.size)) node.size <- 'scaleByExpansion'
  if(missing(max.node.size)) max.node.size <- 20
  if(missing(node.scale.factor)) node.scale.factor <- 6
  if(missing(edge.length)) edge.length <- NULL
  if(missing(edge.width)) edge.width <- NULL
  if(missing(path.list)) path.list <- NULL
  if(missing(specific.node.colors)) specific.node.colors <- NULL
  if(missing(specific.node.shapes)) specific.node.shapes <- NULL
  if(missing(specific.edge.colors)) specific.edge.colors <- NULL
  if(missing(color.by.majority)) color.by.majority <- F
  if(missing(cell.color)) cell.color <- 'seurat_clusters'
  if(missing(specific.cell.colors)) specific.cell.colors <- NULL
  if(missing(cell.size)) cell.size <- 1
  if(missing(network.layout)) network.layout <- 'tree'
  if(missing(save.pdf)) save.pdf <- F
  if(missing(save.dir)) save.dir <- './forests_plots'
  if(missing(show.legend)) show.legend <- T
  if(missing(color.gradient)) color.gradient <- NULL

  bulk.gradient <- NULL

  plot_single_network <- function(g){

    breaks <- 300 #Could try per-clonotype breaks instead of global breaks (or default to ggraph plotting)
    div=2
    # crp <- colorRampPalette(c("#FFFAE8", "#801212"))
    crp <- grDevices::colorRampPalette(colorspace::diverge_hsv(3))
    # length(unique(crp(breaks)))

    #ADD NODE COLORS
    if(is.null(node.color)){

      igraph::V(g)$color[igraph::V(g)$node_type=='sequence'] <- specific.node.colors['sequence']

    }else if(!is.numeric(node.color)){

      features <- igraph::vertex_attr(g, name=node.color, index=which(igraph::V(g)$node_type=='sequence'))

      if(!is.null(igraph::vertex_attr(g, name=paste0(node.color, '_counts'), index=which(igraph::V(g)$node_type=='sequence')))){
        feature_counts <- igraph::vertex_attr(g, name=paste0(node.color, '_counts'), index=which(igraph::V(g)$node_type=='sequence'))

        if(!is.null(color.gradient)) {
          only_bulk_nodes <- c() # finding nodes composed of bulk
          for (i in c(1:length(feature_counts))) {
            if("bulk" %in% names(feature_counts[[i]]) & !("sc" %in% names(feature_counts[[i]])))
              only_bulk_nodes <- append(only_bulk_nodes, i)
          }
          # print(only_bulk_nodes)
        }

        max_indices <- lapply(feature_counts, function(x) which.max(x))
        max_features <- unlist(mapply(function(x,y) x[y], features, max_indices))
      }else{
        max_features <- features
      }

      if(!color.by.majority){
        chosen_colors <- unlist(unname(specific.node.colors[as.character(max_features)]))

        if(!is.null(color.gradient)) {
          for (i in only_bulk_nodes) {

            bulk_count = feature_counts[[i]][["bulk"]]

            if(bulk_count <= 50)
              bulk_count_for_color = bulk_count * 2 # 0 to 100
            else if(bulk_count <= 100)
              bulk_count_for_color = ceiling(50*2/div + bulk_count) # 100 to 150
            else if(bulk_count <= 200)
              bulk_count_for_color = ceiling(50*2/div + 50 + bulk_count/div) # 150 to 200
            else
              bulk_count_for_color = ceiling(50*2/div + 50 + 50*2 - sqrt(5 * 200) + sqrt(5 * bulk_count)) # 200 to ...

            chosen_colors[i] <- crp(breaks)[min(bulk_count_for_color, breaks)]
          }
        }

        igraph::V(g)$color[igraph::V(g)$node_type=='sequence'] <- chosen_colors # in new_plot_t was unlist(chosen_colors)
      }else{
        max_feature_counts <- unlist(lapply(max_features, function(x) length(which(max_features==x))))
        max_features <- max_features[which.max(max_feature_counts)]
        chosen_colors <- unlist(unname(specific.node.colors[as.character(max_features)]))
        igraph::V(g)$color[igraph::V(g)$node_type=='sequence'] <- rep(chosen_colors, length(igraph::V(g)[igraph::V(g)$node_type=='sequence']))
      }

    }else if(is.numeric(node.color)){
      igraph::V(g)$color[igraph::V(g)$node_type=='sequence'] <- '#FFCC00'
      suppressWarnings(igraph::V(g)$color[node.color] <- specific.node.colors[node.color])
    }

    #Temp fix to plot label propagation intermediates

    #### NEXT 2 lines were in new_plot_t, instead of the 19 lines below
    # igraph::V(g)$color[igraph::V(g)$node_type=='intermediate'] <- specific.node.colors['intermediate']
    # igraph::V(g)$color[igraph::V(g)$node_type=='germline'] <- specific.node.colors['germline']

    if(!is.null(node.color)){
      if(!stringr::str_detect(node.color, '_label_propagation')){
        igraph::V(g)$color[igraph::V(g)$node_type=='intermediate'] <- specific.node.colors['intermediate']
        igraph::V(g)$color[igraph::V(g)$node_type=='germline'] <- specific.node.colors['germline']
        if('inferred' %in% names(specific.node.colors)){
          igraph::V(g)$color[igraph::V(g)$node_type=='inferred'] <- specific.node.colors['inferred']
        }
      }else{
        prop_labels <- igraph::vertex_attr(g, name=node.color, index=which(igraph::V(g)$node_type != 'sequence'))
        chosen_colors <- unlist(unname(specific.node.colors[as.character(prop_labels)]))
        igraph::V(g)$color[igraph::V(g)$node_type != 'sequence'] <- chosen_colors
      }
    }else{
      igraph::V(g)$color[igraph::V(g)$node_type=='intermediate'] <- specific.node.colors['intermediate']
      igraph::V(g)$color[igraph::V(g)$node_type=='germline'] <- specific.node.colors['germline']
      if('inferred' %in% names(specific.node.colors)){
        igraph::V(g)$color[igraph::V(g)$node_type=='inferred'] <- specific.node.colors['inferred']
      }
    }


    #ADD NODE SHAPES
    if(is.null(node.shape)){
      igraph::V(g)$shape[igraph::V(g)$node_type=='sequence'] <- 'circle'

    }else if(node.shape=='pie'){
      features <- igraph::vertex_attr(g, name=node.color, index=which(igraph::V(g)$node_type=='sequence'))
      feature_counts <- igraph::vertex_attr(g, name=paste0(node.color, '_counts'), index=which(igraph::V(g)$node_type=='sequence'))
      indices <- lapply(features, function(x) length(x))
      pie_indices <- which(indices!=1)

      for(ind in pie_indices){
        pie_features <- features[[ind]]

        g<-igraph::set_vertex_attr(g, name='shape', index=ind, value='pie')


        if(!is.null(color.gradient)) {
          bulk_count <- feature_counts[ind][[1]][["bulk"]] # we will use it for the color, but we want the node to be 1 half bulk-colored and 1 half sc-colored (see next line)
          feature_counts[ind][[1]][["bulk"]] <- feature_counts[ind][[1]][["sc"]]
        }

        g<-igraph::set_vertex_attr(g, name='pie', index=ind, value=c(feature_counts[ind]))
        # g<-igraph::set_vertex_attr(g, name='pie.color', index=ind, value=list(specific.node.colors[unlist(pie_features)]))

        if(!is.null(color.gradient))
          feature_counts[ind][[1]][["bulk"]] <- bulk_count

        val = list(unlist(specific.node.colors[unlist(features[ind])]))
        # print(val)

        if(!is.null(color.gradient)) {
          if(bulk_count <= 50)
            bulk_count_for_color = bulk_count * 2 # 0 to 100
          else if(bulk_count <= 100)
            bulk_count_for_color = ceiling(50*2/div + bulk_count) # 100 to 150
          else if(bulk_count <= 200)
            bulk_count_for_color = ceiling(50*2/div + 50 + bulk_count/div) # 150 to 200
          else
            bulk_count_for_color = ceiling(50*2/div + 50 + 50*2 - sqrt(5 * 200) + sqrt(5 * bulk_count)) # 200 to ...

          val[[1]][["bulk"]] <- crp(breaks)[min(bulk_count_for_color, breaks)] # ADDED
        }
        # print(feature_counts[ind][[1]][["bulk"]])
        # print(val)
        # print(which(names(val)=="bulk"))
        g<-igraph::set_vertex_attr(g, name='pie.color', index=ind, value=val)
        # val[[1]][["bulk"]] <- "#FF0000" # ADDED, once that we have added the color, we put it back to one so that only a small slide of the pie will have bulk color
        # print(features[ind])
        # print(specific.node.colors[unlist(features[ind])])
        # cat("\n\n")
        # print(names(igraph::V(g)))

        #igraph::V(g)[ind]$shape <- 'pie'
        #igraph::V(g)[ind]$pie <- unlist(feature_counts[ind])
        #igraph::V(g)[ind]$pie.color <- unlist(specific.node.colors[unlist(features[ind])])
      }

      non_pie_indices <- which(indices==1)
      g<-igraph::set_vertex_attr(g, name='shape', index=non_pie_indices, value='circle')

    }else{
      features <- igraph::vertex_attr(g, name=node.shape, index=which(igraph::V(g)$node_type=='sequence'))
      feature_counts <- igraph::vertex_attr(g, name=paste0(node.color, '_counts'), index=which(igraph::V(g)$node_type=='sequence'))
      max_indices <- lapply(feature_counts, function(x) which.max(x))
      max_features <- unlist(mapply(function(x,y) x[y], features, max_indices))
      chosen_shapes <- unlist(unname(specific.node.shapes[max_features]))
      chosen_shapes[is.null(chosen_shapes)] <- 'circle'
      igraph::V(g)$shape[igraph::V(g)$node_type=='sequence'] <- chosen_shapes
    }

    igraph::V(g)$shape[igraph::V(g)$node_type=='intermediate'] <- 'circle'
    igraph::V(g)$shape[is.na(igraph::V(g)$node_type)] <- 'circle'
    igraph::V(g)$shape[igraph::V(g)$node_type=='germline'] <- 'circle'
    #igraph::V(g)$shape[igraph::V(g)$node_type=='cell'] <- 'circle'


    #ADD NODE SIZES
    if(is.null(node.size)){
      igraph::V(g)$size <- node.scale.factor * rep(1, length(igraph::V(g)))
    }else if((node.size %in% c('scaleByExpansion', 'scaleByCloseness', 'scaleByBetweenness'))){
      if(node.size=='scaleByExpansion'){
        #igraph::V(g)$size[igraph::V(g)$size > 100] <- max(igraph::V(g)$size[igraph::V(g)$size > 100]) + 10
        igraph::V(g)$cell_number[is.na(igraph::V(g)$cell_number)] <- 1



        if (!is.null(color.gradient))
          sequence_frequencies <- ifelse(igraph::V(g)$cell_number <7, igraph::V(g)$cell_number, ceiling(sqrt(8 + 5 * igraph::V(g)$cell_number)))  # ADDED, the non-linear function (NEED TO INTRODUCE SOME NON-LINEARITY)
        else
          sequence_frequencies <- ceiling(igraph::V(g)$cell_number/1.5) # was without the division by 1.5
        sequence_frequencies[sequence_frequencies >= max.node.size] <- max.node.size

        #sequence_frequencies[is.na(sequence_frequencies)] <- 1
        scaled_frequencies <- unlist(lapply(sequence_frequencies, function(x) (8+2*x)/10))
        igraph::V(g)$size <- node.scale.factor * scaled_frequencies
      }else if(node.size=='scaleByCloseness'){
        igraph::V(g)$size[igraph::V(g)$node_type!='intermediate'] <- node.scale.factor * igraph::closeness(g, mode = c("all"), weights = NULL, normalized = TRUE)
      }else if(node.size=='scalebyBetweenness'){
        igraph::V(g)$size[igraph::V(g)$node_type!='intermediate'] <- node.scale.factor * igraph::betweenness(g, directed = T, weights = NULL, nobigint = FALSE)
        igraph::V(g)$size[igraph::V(g)$size==0] <- 1
      }else if(node.size=='scaleByEigen'){
        igraph::V(g)$size[igraph::V(g)$node_type!='intermediate'] <- node.scale.factor * igraph::eigen_centrality(g, directed=T, scale=T, weights=NULL)$vector
      }
    }else{
      if(!(node.size %in% igraph::vertex_attr_names(g))){
        stop('Node size feature was not found in your igraph object. Make sure to run VDJ_networks with node.features=node size feature')
      }
      igraph::V(g)$size[igraph::V(g)$node_type=='sequence'] <- node.scale.factor * igraph::vertex_attr(g, name=node.size, index=which(igraph::V(g)$node_type=='sequence'))
    }

    igraph::V(g)$size[igraph::V(g)$node_type=='germline'] <- min(igraph::V(g)$size[igraph::V(g)$node_type=='sequence']) * 1.5 # was 0.8
    igraph::V(g)$size[igraph::V(g)$node_type=='intermediate'] <- min(igraph::V(g)$size[igraph::V(g)$node_type=='sequence']) * 0.8


    if(is.null(node.label)){
      igraph::V(g)$label <- NA

    }else if(node.label == 'cells'){
      igraph::V(g)$cell_number[is.na(igraph::V(g)$cell_number)] <- 1
      igraph::V(g)$label <- ifelse(igraph::V(g)$cell_number>1, igraph::V(g)$cell_number, "") # We don't want to show all the 1s
      if (!is.null(color.gradient))
        igraph::V(g)$label <- as.numeric(igraph::V(g)$label) - 1 # ADDED: FOR BULK AND SC MERGED TREES ONLY
      igraph::V(g)$label[igraph::V(g)$node_type == 'germline'] <- NA
      igraph::V(g)$label[igraph::V(g)$node_type == 'intermediate'] <- NA
      igraph::V(g)$label[igraph::V(g)$node_type == 'bulk'] <- NA
      igraph::V(g)$label[igraph::V(g)$node_type == 'inferred'] <- NA


    }else{
      igraph::V(g)$label <- 1:length(igraph::V(g))
    }

    label_size <- rep(min(igraph::V(g)$size)*0.150, length(igraph::V(g))) # was 0.125
    igraph::V(g)$label.cex <- label_size
    igraph::V(g)$label.font <- 2

    if (bulk.gradient == T)
      igraph::V(g)$label.color <- '#F9973B' #'#8FA5A6' #'#FFFFCC'
    else
      igraph::V(g)$label.color <- 'beige'

    #Quick fix
    igraph::V(g)$shape[is.na(igraph::V(g)$shape)] <- 'circle'


    #ADD EDGE LENGTHS
    if(!is.null(edge.length)){
      if(edge.length=='scaleByWeight'){
        igraph::E(g)$length <- 1/igraph::E(g)$weight
      }else{
        igraph::E(g)$length <- igraph::edge_attr(g, name=edge.length)
      }
    }

    #ADD EDGE WIDTHS
    if(!is.null(edge.width)){
      if(edge.width=='scaleByWeight'){
        igraph::E(g)$width <- 1/igraph::E(g)$weight
      }else{
        igraph::E(g)$width <- igraph::edge_attr(g, name=edge.width)
      }
    }

    if(graph.type!= 'heterogeneous'){
      igraph::E(g)$weight <- 1
    }

    #COLOR PATHS
    #igraph::E(g)$color <- 'grey24'
    #if(!is.null(path.list)){
    #  path_names <- names(paths)
    #  if(length(paths)==0){
    #    igraph::E(g)$color <- 'grey24'
    #  }else{
    #    for(i in 1:length(paths)){
    #      if(length(paths[[i]])==0){
    #        igraph::E(g)$color <- 'grey24'
    #      }else{
    #        for(j in 1:length(paths[[i]])){
    #          print(g)
    #          print(paths[[i]][[j]])
    #          igraph::E(g, path=paths[[i]][[j]])$color <- specific.edge.colors[path_names[i]]
    #        }
    #      }
    #    }
    #  }
    #}

    igraph::E(g)$arrow.size <- min(igraph::V(g)$size) * 0.065


    if(network.layout=='tree'){
      layout <- igraph::layout_as_tree(g, root=which(igraph::V(g)$node_type=='germline'), circular=F)
    }else if(network.layout=='fr'){
      layout <- igraph::layout_with_fr(g)
    }else{
      layout <- igraph::layout_nicely(g)
    }

    g$layout <- layout

    return(g)
  }

  plot_single_network_cells <- function(g){

    g <- igraph::set_vertex_attr(g, 'lvl', index = igraph::V(g), ifelse(igraph::V(g)$type == 'sequence', 1, 2))

    igraph::V(g)$shape[igraph::V(g)$node_type=='cell'] <- 'square'
    igraph::V(g)$size[igraph::V(g)$node_type=='cell'] <- node.scale.factor * cell.size
    igraph::V(g)$label[igraph::V(g)$node_type=='cell'] <- NA


    if(is.null(cell.color)){
      igraph::V(g)$color[igraph::V(g)$node_type=='cell'] <- unname(unlist(specific.cell.colors))
    }else{
      igraph::V(g)$color[igraph::V(g)$node_type=='cell'] <- specific.cell.colors[unlist(igraph::vertex_attr(g, name = cell.color, index = which(igraph::V(g)$node_type=='cell')))]
    }


    igraph::E(g)$width[igraph::E(g)$type == 'cell'] <- 0.8 * igraph::E(g)$weight[igraph::E(g)$type == 'cell']
    igraph::E(g)$width[igraph::E(g)$type == 'interlevel'] <- 0.4 * igraph::E(g)$weight[igraph::E(g)$type == 'interlevel']
    igraph::E(g)$width[igraph::E(g)$type == 'sequence'] <- 1.5

    igraph::E(g)$color[igraph::E(g)$type == 'cell'] <- 'gray'
    igraph::E(g)$color[igraph::E(g)$type == 'sequence'] <- 'black'
    igraph::E(g)$color[igraph::E(g)$type == 'interlevel'] <- 'gray'

    igraph::E(g)$arrow.size[igraph::E(g)$type == 'sequence'] <- min(igraph::V(g)$size[igraph::V(g)$node_type == 'sequence']) * 0.065
    #Quick fix for pie shapes in heteogeneous expanded
    igraph::V(g)$shape[is.na(igraph::V(g)$shape)] <- 'circle'


    #igraph::E(g)$color[igraph::E(g)$type == 'sequence'] <- 'black'




    if(network.layout=='tree'){
      layout <- purrr::partial(igraph::layout_as_tree, root = which(igraph::V(g)$node_type == 'germline'))
    }else if(network.layout=='fr'){
      layout <- purrr::partial(igraph::layout_with_fr)
    }else{
      layout <- purrr::partial(igraph::layout_nicely)
    }


    xy <- graphlayouts::layout_as_multilevel(g,
                                             type = "separate",
                                             FUN1 = layout,
                                             FUN2 = graphlayouts::layout_with_stress,
                                             alpha = 25, beta = 45
    )

    g$layout <- xy

    return(g)
  }


  standard_colors <- c('germline' = 'white', 'intermediate' = 'grey86')
  sequence_colors <- c('sequence' = '#FFCC00')
  cell_colors <- c('cell' = 'bisque2')

  shapes <- setdiff(igraph::shapes(), "")
  shapes <- shapes[shapes!='pie' & shapes!='none']

  if(!is.null(node.color)){
    unique_features <- c()

    if(inherits(network.list, 'list')){
      for(i in 1:length(network.list)){
        for(j in 1:length(network.list[[i]])){
          g <- network.list[[i]][[j]]@tree
          unique_features <- unlist(unique(c(unique_features, igraph::vertex_attr(g, name=node.color, index=which(igraph::V(g)$node_type=='sequence')))))
        }
      }
    }else{
      g <- network.list@tree
      unique_features <- unlist(unique(igraph::vertex_attr(g, name=node.color, index=which(igraph::V(g)$node_type=='sequence'))))
    }

    unique_features <- unique(unique_features)

    if(length(unique_features) == 0){
      message('Found 0 unique features for node colors, will use the default node colors instead!')
      specific.node.colors <- c(sequence_colors, standard_colors)
    }

    if(is.null(specific.node.colors)){
      if(!is.numeric(unique_features)){
        unique_feature_colors <- c(grDevices::rainbow(length(unique_features)))

      }else{
        breaks <- length(unique_features) #Could try per-clonotype breaks instead of global breaks (or default to ggraph plotting)
        pal <- grDevices::colorRampPalette(c('cornflowerblue','darkblue'))
        unique_feature_colors <- pal(breaks)[as.numeric(cut(unique_features, breaks = breaks))]
      }

      names(unique_feature_colors) <- unique_features
      specific.node.colors <- c(unique_feature_colors, standard_colors)
      specific.node.colors <- specific.node.colors[order(nchar(names(specific.node.colors)),names(specific.node.colors))]


    }else if(!inherits(specific.node.colors, 'list')){
      if(!is.numeric(unique_features)){
        features_not_added <- setdiff(unique_features, names(specific.node.colors))
        if(length(features_not_added) > 0){
          colors_to_add <- c(grDevices::rainbow(length(features_not_added)))
          message(paste0('Colors not added for the following unique features in the specific.node.colors parameter:  ', paste0(features_not_added, collapse = ','), '.Will default to the following colors: ', paste0(colors_to_add, collapse = ',')))
          names(colors_to_add) <- features_not_added
          specific.node.colors <- c(specific.node.colors, colors_to_add)
        }
        specific.node.colors <- c(specific.node.colors, standard_colors)

      }else{
        if(length(specific.node.colors) == 1){
          c1 <- 'gray'
          c2 <- unlist(specific.node.colors)
        }else{
          c1 <- specific.node.colors[1]
          c2 <- specific.node.colors[2]
        }

        breaks <- length(unique_features) #Could try per-clonotype breaks instead of global breaks (or default to ggraph plotting)
        pal <- grDevices::colorRampPalette(c(c1, c2))
        unique_feature_colors <- pal(breaks)[as.numeric(cut(unique_features, breaks = breaks))]
        names(unique_feature_colors) <- unique_features
        specific.node.colors <- unique_feature_colors
      }
    }else if(inherits(specific.node.colors, 'list')){
      if(!('germline' %in% names(specific.node.colors))){
        germ <- list('blue')
        names(germ) <- 'germline'

        specific.node.colors <- c(specific.node.colors, germ)
      }

      if(!('intermediate' %in% names(specific.node.colors))){
        interm <- list('grey86')
        names(interm) <- 'intermediate'

        specific.node.colors <- c(specific.node.colors, germ)

      }
      specific.node.colors <- specific.node.colors[order(nchar(names(specific.node.colors)),names(specific.node.colors))]

    }

  }else{
    specific.node.colors <- c(specific.node.colors, sequence_colors, standard_colors)
    specific.node.colors <- specific.node.colors[order(nchar(names(specific.node.colors)),names(specific.node.colors))]
  }

  if('unknown' %in% names(specific.node.colors)){
    specific.node.colors['unknown'] <- 'burlywood'
  }



  if(is.null(specific.node.shapes) & !is.null(node.shape)){
    if(node.shape!='pie'){
      unique_features <- NULL
      if(inherits(network.list, 'list')){
        for(i in 1:length(network.list)){
          for(j in 1:length(network.list[[i]])){
            g <- network.list[[i]][[j]]@tree
            unique_features <- unlist(unique(c(unique_features, igraph::vertex_attr(g, name=node.shape, index=which(igraph::V(g)$node_type=='sequence')))))
          }
        }
      }else{
        g <- network.list@tree
        unique_features <- unlist(unique(c(unique_features, igraph::vertex_attr(g, name=node.shape, index=which(igraph::V(g)$node_type=='sequence')))))
      }

      unique_features <- unique(unique_features)

      if(length(shapes) < length(unique_features)){
        stop('Not enough shapes for the number of unique features in the networks')
      }

      specific.node.shapes <- as.list(shapes[1:length(unique_features)])
      names(specific.node.shapes) <- unique_features
    }
  }

  #if(!is.null(path.list)){
  #  if(is.null(specific.edge.colors)){
  #    path_names <- NULL
  #    for(i in 1:length(path.list)){
  #      for(j in 1:length(path.list[[i]])){
  #        path_names <- unique(c(path_names,names(path.list[[i]][[j]])))
  #      }
  #    }
  #    specific.edge.colors <- c(grDevices::rainbow(length(path_names), rev=T))
  #    names(specific.edge.colors) <- path_names
  #  }
  #}

  if(graph.type == 'heterogeneous'){

    if(!is.null(cell.color)){
      if(inherits(network.list, 'list')){
        unique_features <- c()
        for(i in 1:length(network.list)){
          for(j in 1:length(network.list[[i]])){
            g <- network.list[[i]][[j]]@heterogeneous
            unique_features <- unlist(unique(c(unique_features, igraph::vertex_attr(g, name=cell.color, index=which(igraph::V(g)$node_type=='cell')))))
          }
        }
      }else{
        g <- network.list@heterogeneous
        unique_features <- unlist(unique(igraph::vertex_attr(g, name=cell.color, index=which(igraph::V(g)$node_type=='cell'))))

      }

      if(cell.color == 'seurat_clusters'){
        unique_features <- unlist(lapply(unique_features, function(x) as.factor(x)))
      }

      unique_features <- unique(unique_features)

      if(length(unique_features) == 0){
        message('Found 0 unique features for cell colors, will use the default node colors instead!')
        specific.cell.colors <- cell_colors

      }else{
        if(is.null(specific.cell.colors)){
          if(!is.numeric(unique_features)){
            unique_feature_colors <- c(grDevices::topo.colors(length(unique_features)))

          }else{
            breaks <- length(unique_features) #Could try per-clonotype breaks instead of global breaks (or default to ggraph plotting)
            pal <- grDevices::colorRampPalette(c('darkolivegreen1','darkolivegreen4'))
            unique_feature_colors <- pal(breaks)[as.numeric(cut(unique_features, breaks = breaks))]
          }

          names(unique_feature_colors) <- unique_features
          specific.cell.colors <- unique_feature_colors
          specific.cell.colors <- specific.cell.colors[order(names(specific.cell.colors), nchar(names(specific.cell.colors)))]

        }else{
          if(!is.numeric(unique_features)){
            features_not_added <- setdiff(unique_features, names(specific.cell.colors))
            if(length(features_not_added) > 0){
              colors_to_add <- c(grDevices::topo.colors(length(features_not_added)))
              message(paste0('Colors not added for the following unique features in the specific.cell.colors parameter:  ', paste0(features_not_added, collapse = ','), '.Will default to the following colors: ', paste0(colors_to_add, collapse = ',')))
              names(colors_to_add) <- features_not_added
              specific.cell.colors <- c(specific.cell.colors, colors_to_add)
            }

          }else{
            if(length(specific.cell.colors) == 1){
              c1 <- 'gray'
              c2 <- unlist(specific.cell.colors)
            }else{
              c1 <- specific.cell.colors[1]
              c2 <- specific.cell.colors[2]
            }

            breaks <- length(unique_features) #Could try per-clonotype breaks instead of global breaks (or default to ggraph plotting)
            pal <- grDevices::colorRampPalette(c(c1, c2))
            unique_feature_colors <- pal(breaks)[as.numeric(cut(unique_features, breaks = breaks))]
            names(unique_feature_colors) <- unique_features
            specific.cell.colors <- unique_feature_colors
          }
        }
      }
    }else{
      specific.cell.colors <- cell_colors
    }
  }

  if(inherits(network.list, 'list')){

    for(i in 1:length(network.list)){

      for(j in 1:length(network.list[[i]])){

        if(graph.type == 'tree'){
          network.list[[i]][[j]]@plot_ready <- network.list[[i]][[j]]@tree %>%
            plot_single_network()

        }else if(graph.type == 'heterogeneous'){
          network.list[[i]][[j]]@plot_ready <- network.list[[i]][[j]]@heterogeneous %>%
            plot_single_network() %>%
            plot_single_network_cells()


        }else if(graph.type == 'dynamic'){
          network.list[[i]][[j]]@plot_ready <- network.list[[i]][[j]]@dynamic %>%
            plot_single_network()



        }else{
          stop('Unrecognized graph type!')
        }

        title_name <- paste0(network.list[[i]][[j]]@sample_id, '_', network.list[[i]][[j]]@clonotype_id)

        cex_ = 0.75
        inset_ = -0.14


        if(save.pdf){
          if(!dir.exists(save.dir)) dir.create(save.dir)
          file_name <- paste0(title_name, '.pdf')

          grDevices::pdf(paste0(save.dir, '/', file_name))
          plot(network.list[[i]][[j]]@plot_ready)
          graphics::title(title_name)
          graphics::legend(
            "left",
            legend = unlist(names(specific.node.colors)),
            pt.bg  = unlist(specific.node.colors),
            pch    = 21,
            cex    = cex_,
            inset = inset_,
            bty    = "n",
            title  = node.color
          )

          if(!is.null(specific.cell.colors)){
            graphics::legend(
              "right",
              legend = unlist(names(specific.cell.colors)),
              pt.bg  = unlist(specific.cell.colors),
              pch    = 21,
              cex    = 0.5,
              bty    = "n",
              title  = cell.color
            )
          }
          grDevices::dev.off()
        }else{

          plot(network.list[[i]][[j]]@plot_ready)
          graphics::title(title_name)
          if(show.legend){
            graphics::legend(
              "left",
              legend = unlist(names(specific.node.colors)),
              pt.bg  = unlist(specific.node.colors),
              pch    = 21,
              cex    = 0.5,
              bty    = "n",
              title  = node.color
            )

            if(!is.null(specific.cell.colors)) {
              graphics::legend(
                "right",
                legend = unlist(names(specific.cell.colors)),
                pt.bg  = unlist(specific.cell.colors),
                pch    = 21,
                cex    = 0.5,
                bty    = "n",
                title  = cell.color
              )
            }
          }

        }
      }
    }

  }else{

    if(graph.type == 'tree'){
      network.list@plot_ready <- network.list@tree %>%
        plot_single_network()

    }else if(graph.type == 'heterogeneous'){
      network.list@plot_ready <- network.list@heterogeneous %>%
        plot_single_network() %>%
        plot_single_network_cells()


    }else if(graph.type == 'dynamic'){
      network.list@plot_ready <- network.list@dynamic %>%
        plot_single_network()


    }else{
      stop('Unrecognized graph type!')
    }


    title_name <- paste0(network.list@sample_id, '_', network.list@clonotype_id)
    if(save.pdf){
      if(!dir.exists(save.dir)) dir.create(save.dir)
      file_name <- paste0(title_name, '.pdf')

      grDevices::pdf(paste0(save.dir, '/', file_name))

      plot(network.list@plot_ready)
      #title(title_name)
      if(show.legend){
        graphics::legend(
          "left",
          legend = unlist(names(specific.node.colors)),
          pt.bg  = unlist(specific.node.colors),
          pch    = 21,
          cex    = 0.5,
          bty    = "n",
          title  = node.color
        )


        if(!is.null(specific.cell.colors)){
          graphics::legend(
            "right",
            legend = unlist(names(specific.cell.colors)),
            pt.bg  = unlist(specific.cell.colors),
            pch    = 21,
            cex    = 0.5,
            bty    = "n",
            title  = cell.color
          )
        }
      }

      grDevices::dev.off()
    }else{

      plot(network.list@plot_ready)
      #title(title_name)

      if(show.legend){
        graphics::legend(
          "left",
          legend = unlist(names(specific.node.colors)),
          pt.bg  = unlist(specific.node.colors),
          pch    = 21,
          cex    = 0.5,
          bty    = "n",
          title  = node.color
        )


        if(!is.null(specific.cell.colors)){
          graphics::legend(
            "right",
            legend = unlist(names(specific.cell.colors)),
            pt.bg  = unlist(specific.cell.colors),
            pch    = 21,
            cex    = 0.5,
            bty    = "n",
            title  = cell.color
          )
        }
      }
    }
  }

  return(network.list)
}
