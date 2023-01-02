#' Node degree bar plots colored by sequence-level features


#'@description Creates an ordered bar plot of the vertices with the highest degree, colored by the sequence-level features specified in the features parameter. Creates a bar plot per sample, per clonotype, or per entire AntibodyForests object, depending on the input format (object, list, nested list).
#' Secondary option to create a node degree histogram if plot.type == 'histogram'.
#' @param trees AntibodyForests object/list of AntibodyForests objects - the resulting sequence similarity or minimum spanning tree networks from the AntibodyForests function
#' @param features vector of strings - sequence-level features for the bar plot colors.
#' @param graph.type string - the graph type available in the AntibodyForests object which will be used as the function input.
#' Currently supported network/analysis types: 'tree' (for the minimum spanning trees or sequence similarity networks obtained from the main AntibodyForests function), 'heterogeneous' for the bipartite graphs obtained via AntibodyForests_heterogeneous.
#' @param which.bipartite string - whether to create a bar plot for the cell layer of the bipartite/heterogeneous graph ('cells'), sequence layer ('sequences') or for both ('both').
#' @param max.nodes integer - the maximum number of nodes (bars) to be plotted in the degree bar plot
#' @param nbins integer - number of bins for the degree histograms if plot.type == 'histogram'
#' @param exclude.germline boolean - if T, will exclude the germline nodes when creating the bar plots.
#' @param exclude.intermediates boolean - if T, will exclude the intermediate nodes when creating the bar plots (if expand.intermediates was set to T when creating the AntibodyForests object).
#' @param parallel boolean - if T, will execute the main loop in parallel using mclapply.
#' @param plot.type string - 'barplot' for ordered node degree bar plots colored by the sequence-level feature specified in the features parameter. 'histogram' for node degree histograms.

#' @return a ggplot2 bar plot or histogram.
#' @export
#' @seealso AntibodyForests, AntibodyForests_metrics
#' @examples
#' \dontrun{
#' AntibodyForests_object |> AntibodyForests_degrees(features = 'VDJ_cgene', plot.type = 'barplot')
#'}


AntibodyForests_degrees <- function(trees,
                                    features,
                                    graph.type,
                                    which.bipartite,
                                    max.nodes,
                                    nbins,
                                    exclude.germline,
                                    exclude.intermediates,
                                    parallel,
                                    plot.type){

  if(missing(trees)) stop('Please input an AntibodyForests object or a list of objects!')
  if(missing(features)) features <- NULL
  if(missing(graph.type)) graph.type <- 'tree'
  if(missing(which.bipartite)) which.bipartite <- 'combined'
  if(missing(max.nodes)) max.nodes <- 50
  if(missing(nbins)) nbins <- 20
  if(missing(exclude.germline)) exclude.germline <- T
  if(missing(exclude.intermediates)) exclude.intermediates <-  T
  if(missing(parallel)) parallel <- F
  if(missing(plot.type)) plot.type <- 'barplot' #Or histogram


  get_feature_names <- function(trees, features){

   if(is.null(features)){
     if(inherits(trees, 'list')){
       features <- trees[[1]][[1]]@feature_names
     }else if(inherits(trees, 'AntibodyForests')){
       features <- trees@feature_names
     }

     if(is.null(features)){
       stop('Could not find the features to perform label propagation on! Please provide the feature names in the features parameter!')
     }
   }

   return(features)
  }

  get_graph <- function(tree){

    if(graph.type == 'tree'){
      g <- list(tree@tree)
      names(g) <- 'sequence'

    }else if(graph.type == 'heterogeneous'){
      g <- tree@heterogeneous

      cell_vertices <- which(igraph::V(g)$type == 'cell')
      sequence_vertices <- which(igraph::V(g)$type == 'sequence')
      sequence_g <- igraph::delete_vertices(g, cell_vertices)
      cell_g <- igraph::delete_vertices(g, sequence_vertices)

      if(which.bipartite == 'sequences'){
        g <- list(sequence_g)
      }else if(which.bipartite == 'cells'){
        g <- list(cell_g)
      }else if(which.bipartite == 'combined'){
        g <- list(g)
      }else{
        g <- list()
        g[[1]] <- sequence_g
        g[[2]] <- cell_g
        names(g) <- c('sequence', 'cell')
      }

    }else if(graph.type == 'dynamic'){
      g <- list(tree@dynamic)
      names(g) <- 'sequence'

    }else{
      stop('Graph type not found!')
    }

    if(is.null(g)){
      stop(paste0('Could not find the ', graph.type, ' graph!'))
    }

    return(g)
  }

  get_degrees <- function(g, trees, feature_names){

    for(i in 1:length(g)){
      g[[i]] <- igraph::set_vertex_attr(graph = g[[i]], name = 'degree' , value = igraph::degree(g[[i]]))
    }

    if(length(g) == 2){
      g <- assemble_bipartite(g, original_graph)
    }else{
      g <- g[[1]]
    }

    node_df <- igraph::as_data_frame(g, what = 'vertices')

    for(col in feature_names){
      if(paste0(col, '_counts') %in% colnames(node_df)){
        max_indices <- node_df[paste0(col, '_counts')][which(node_df$node_type=='sequence'),]
        max_indices <- unlist(lapply(max_indices, function(x) which.max(x)))
        max_values <- unlist(mapply(function(x,y) x[y], node_df[col][which(node_df$node_type=='sequence'),], max_indices))
        node_df[col][which(node_df$node_type=='sequence'),] <- max_values
        node_df[col] <- unlist(node_df[col]) #to also change the column class from list into numeric/character

        node_df <- node_df[,!(names(node_df) %in% paste0(col, '_counts'))]

      }
    }

    if(exclude.germline){
      node_df <- node_df[node_df$node_type != 'germline',]
    }

    if(exclude.intermediates){
      node_df <- node_df[node_df$node_type != 'intermediate',]
    }

    degrees <- unlist(node_df$degree)
    node_df <- node_df[order(degrees, decreasing = T),]
    node_df$rank <- 1:nrow(node_df)

    if(!is.null(max.nodes)){
      if(max.nodes < nrow(node_df)){
        node_df <- node_df[1:max.nodes,]
      }
    }

    return(node_df)
  }

  assemble_bipartite <- function(graphs, trees){
    original_graph <- trees@heterogeneous
    sequence_g <- graphs[[1]]
    cell_g <- graphs[[2]]

    sequence_vertices <- igraph::as_data_frame(sequence_g, what = 'vertices')
    cell_vertices <- igraph::as_data_frame(cell_g, what = 'vertices')

    edgelist <- igraph::as_edgelist(original_graph)
    g <- igraph::graph_from_data_frame(d = edgelist, directed = F, vertices = rbind(sequence_vertices, cell_vertices))

    return(g)

  }

  degree_barplot <- function(degree_df, feature_names){
    out_plots <- list()
    for(i in 1:length(feature_names)){
      feature <- feature_names[i]

      out_plots[[i]] <- ggplot2::ggplot(degree_df, ggplot2::aes_string(fill = feature, y = 'degree', x = 'rank')) +
                   ggplot2::geom_bar(stat="identity", width=1, color="black") +
                   ggplot2::scale_y_continuous(expand = c(0,0)) +
                   ggplot2::theme_bw() +
                   ggplot2::theme_classic() +
                   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x =  ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
                   ggplot2::labs(title = feature_names[i]) +
                   ggplot2::labs(x = "Degree rank", y = "Degree")
    }
    names(out_plots) <- feature_names

    return(out_plots)
  }

  degree_histogram <- function(degree_df){

    out_plot <- ggplot2::ggplot(degree_df, ggplot2::aes(x = degree)) +
                   ggplot2::geom_histogram(ggplot2::aes(y=..density..), colour="black", fill="lightblue", bins = nbins) +
                   ggplot2::theme_bw() +
                   ggplot2::theme_classic() +
                   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x =  ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
                   ggplot2::labs(title = 'Degree histogram') +
                   ggplot2::labs(x = "Degree", y = "Density")

    return(out_plot)
  }

  get_plot <- function(tree, features){
    g <- get_graph(tree)
    degree_df <- get_degrees(g, tree, features)

    if(plot.type == 'histogram'){
      p <- degree_histogram(degree_df)

    }else{
      p <- degree_barplot(degree_df, features)
    }

    return(p)
  }

  features <- get_feature_names(trees, features)
  if(inherits(trees, 'list')){
    plots <- list()
    for(i in 1:length(trees)){

      if(parallel){
        requireNamespace('parallel')
        cores <- parallel::detectCores()
        plots[[i]] <- parallel::mclapply(trees[[i]], mc.cores = cores,
                                                 FUN = function(x) {get_plot(x, features)
                                                                    })

      }else{
        plots[[i]] <- lapply(trees[[i]], function(x) get_plot(x, features))

      }
      names(plots[[i]]) <- names(trees[[i]])
    }
    names(plots) <- names(trees)

  }else if(inherits(trees, 'AntibodyForests')){
    plots <- get_plot(trees, features)

  }else{
    stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
  }


  return(plots)
}
