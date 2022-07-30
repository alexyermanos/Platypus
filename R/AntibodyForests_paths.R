#' Calculates the longest/shortest paths from a node to a given node for the AntibodyForests minimum spanning trees / sequence similarity networks

#'@description Calculates the longest or shortest paths in a given AntibodyForests graph from a node to another given node. Nodes can be specified as integers (e.g., path.from = 5, picking the fifth node in the igraph vertex list) or by predetermined attributes (e.g., path.from = 'germline' and path.to = 'hub' will calculate the paths between all germlines and hubs in a set of networks.
#' Moreover, there is an option to select paths for nodes given specific node features (e.g., path.from = list('seurat_clusters', '1') and path.to = 'hub' will infer the paths from the nodes with a majority of Seurat clusters = 1, and to the hub nodes).

#' @param trees nested list of AntibodyForests objects or single object, as obtained from the AntibodyForests function.
#' @param graph.type string - the graph type available in the AntibodyForests object which will be used as the function input.
#' Currently supported network/analysis types: 'tree' (for the minimum spanning trees or sequence similarity networks obtained from the main AntibodyForests function), 'heterogeneous' for the bipartite graphs obtained via AntibodyForests_heterogeneous, 'dynamic' for the dynamic networks obtained from AntibodyForests_dynamics.
#' @param path.from string/list of strings/integer - starting nodes for a path. Options are either an integer, selecting the nth most abundant node, 'hub' to select the hub nodes, 'most_expanded' for the nodes with most cells, 'leaf' for leaf nodes, 'germline' for the germline node, or a list of the form list(feature, feature_value) to select nodes with a specific feature value.
#' @param path.to string/list of strings/integer - end nodes for a path. Options are either an integer, selecting the nth most abundant node, 'hub' to select the hub nodes, 'most_expanded' for the nodes with most cells, 'leaf' for leaf nodes, 'germline' for the germline node, or a list of the form list(feature, feature_value) to select nodes with a specific feature value.
#' @param paths string - whether to calculate the longest path ('longest'), shortest path ('shortest'), or both (c('longest', 'shortest'))
#' @param interlevel.from  string/list of strings/integer - starting nodes for an interlevel path for the bipartite/heterogeneous networks (if graph.type = 'heterogeneous'). Options are either an integer, selecting the nth most abundant node, 'hub' to select the hub nodes, 'most_expanded' for the nodes with most cells, 'leaf' for leaf nodes, 'germline' for the germline node, or a list of the form list(feature, feature_value) to select nodes with a specific feature value.
#' @param weighted boolean - whether to calculate the weighted or unweighted shortest/longest paths.
#' @param plot.results boolean - if T, will output a bar plot of node feature counts per path (of all nodes in a given path). Features are determined by the color.by parameter.
#' @param color.by string - features for the feature count per path bar plots if plot.results is set to T.
#' @param cell.frequency boolean - whether to consider the node or cell frequency for the feature counts in the resulting bar plot if plot.results is T.
#' @param parallel boolean - whether to execute the main subroutine in parallel or not. Requires the 'parallel' R package to be installed.


#' @return nested list of AntibodyForests objects or a single object with a new slot - paths. If plot.results is T, will also output a bar plot of feature counts per path (considering all node in a given path).
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_paths(trees, graph.type = 'tree',
#' path.from = 'germline', path.to = 'leaf',
#' plot.results = T, color.by = 'seurat_clusters')
#'}


AntibodyForests_paths <- function(trees,
                                  graph.type,
                                  path.from,
                                  path.to,
                                  paths,
                                  interlevel.from,
                                  weighted,
                                  plot.results,
                                  color.by,
                                  cell.frequency,
                                  parallel
                                  ){

  if(missing(trees)) stop('Please input a nested list of AntibodyForests objects to get the paths for!')
  if(missing(graph.type)) graph.type <- 'tree'
  if(missing(path.from)) path.from <- 'germline'
  if(missing(path.to)) path.to <- 'leaf'
  if(missing(paths)) paths <- c('shortest')
  if(missing(interlevel.from)) interlevel.from <- 'hub'
  if(missing(weighted)) weighted <- T
  if(missing(plot.results)) plot.results <- F
  if(missing(color.by)) color.by <- 'VDJ_cgene'
  if(missing(cell.frequency)) cell.frequency <- F
  if(missing(parallel)) parallel <- T

  features <- NULL
  counts <- NULL
  clonotype_id <- NULL
  total_cells <- NULL
  total_nodes <- NULL
  sample_id <- NULL

  get_paths <- function(tree){

    if(graph.type == 'tree'){
      g <- tree@tree

    }else if(graph.type == 'heterogeneous'){
      het_g <- tree@heterogeneous
      if(is.null(het_g)){
        stop('Please call AntibodyForests_heterogeneous before for interlevel graph paths!')
      }

      cell_nodes <- which(igraph::V(het_g)$node_type == 'cell')
      g <- igraph::delete_vertices(het_g, cell_nodes)


    }

    sequence_nodes <- which(igraph::V(g)$node_type=='sequence')
    germline_node <-  which(igraph::V(g)$node_type=='germline')

    node_list <- list(path.from, path.to, interlevel.from)
    final_nodes <- list()
    output_paths <- c()

    for(i in 1:length(node_list)){
      if(!is.null(node_list[[i]])){
        if(is.list(node_list[[i]])){
          feature <- unlist(node_list[[i]])[1]
          value <- unlist(node_list[[i]])[2]
          max_indices <- igraph::vertex_attr(g, name=paste0(feature, '_counts'), index=which(igraph::V(g)$node_type=='sequence'))
          max_indices <- lapply(max_indices, function(x) which.max(x))
          selected_max_feature_values <- unlist(mapply(function(x,y) x[y], igraph::vertex_attr(g, name=feature, index=which(igraph::V(g)$node_type=='sequence')), max_indices))
          g<-igraph::set_vertex_attr(g, name=feature, index=which(igraph::V(g)$node_type=='sequence'), value=selected_max_feature_values)

          if(is.numeric(unlist(node_list[[i]])[2])){
            final_nodes[[i]] <- which(igraph::vertex_attr(g, name=feature, index=sequence_nodes) > value)

          }else{
            final_nodes[[i]] <- which(igraph::vertex_attr(g, name=feature, index=sequence_nodes) %in% value)
          }

          if(length(node_list[[i]])==3){
            if(unlist(node_list[[i]])[3]=='leaf'){
              final_nodes[[i]] <- base::intersect(unlist(final_nodes[[i]]), which(igraph::V(g)$leaf=='leaf'))
            }else if(unlist(node_list[[i]])[3]=='most_expanded'){
              final_nodes[[i]] <- base::intersect(unlist(final_nodes[[i]]), which(igraph::V(g)$most_expanded=='yes'))
            }else if(unlist(node_list[[i]][3]=='hub')){
              final_nodes[[i]] <- base::intersect(unlist(final_nodes[[i]]), which(igraph::V(g)$hub=='yes'))
            }
          }

        }else if(is.numeric(node_list[[i]])){
          final_nodes[[i]] <- node_list[[i]]
        }else if(node_list[[i]]=='all'){
          final_nodes[[i]] <- c(sequence_nodes, germline_node)
        }else if(node_list[[i]]=='leaf'){
          final_nodes[[i]] <- which(igraph::V(g)$leaf=='leaf')
        }else if(node_list[[i]]=='germline'){
          final_nodes[[i]] <- germline_node
        }else if(node_list[[i]]=='most_expanded'){
          final_nodes[[i]] <- which(igraph::V(g)$most_expanded=='yes')
        }else if(node_list[[i]]=='hub'){
          final_nodes[[i]] <- which(igraph::V(g)$hub=='yes')
        }else{
          return(NULL)
        }
      }else{
        return(NULL)
      }
    }


    from_nodes <- unlist(final_nodes[[1]])
    to_nodes <- unlist(final_nodes[[2]])
    interlevel_from <- unlist(final_nodes[[3]])


    if(length(from_nodes)==0 | length(to_nodes)==0){
      return(NULL)
    }

    output_paths <- c()
    for(i in 1:length(paths)){
      path_type <- paths[i]
      for(node in from_nodes){
        #to_nodes <- to_nodes[to_nodes!=node]
        if(path_type != 'interlevel'){
          if(path_type!='all'){
            if(weighted){
              distances <- suppressWarnings(igraph::distances(g, v=node, to=to_nodes, weights=NULL))
              #non_inf_ind <- which(distances != 'Inf')
              #final_to_nodes <- as.integer(colnames(distances))[non_inf_ind]
              #distances <- distances[non_inf_ind]
              final_to_nodes <- to_nodes

            }else{
              distances <- suppressWarnings(igraph::distances(g, v=node, to=to_nodes, weights=NA))
              #non_inf_ind <- which(distances != 'Inf')
              #final_to_nodes <- as.integer(colnames(distances))[non_inf_ind]
              #distances <- distances[non_inf_ind]
              final_to_nodes <- to_nodes
            }

            if(path_type=='longest'){
              final_to_nodes <- final_to_nodes[which(distances==max(distances))]
            }else if(path_type=='shortest'){
              final_to_nodes <- final_to_nodes[which(distances==min(distances))]

            }else{
              stop(paste0('Incorrect path argument - ', path_type))
            }
           }

            if(weighted){
              out <- suppressWarnings(igraph::all_shortest_paths(g, from=node, to=final_to_nodes, weights=NULL))
            }else{
              out <- suppressWarnings(igraph::all_shortest_paths(g, from=node, to=final_to_nodes, weights=NA))
            }
            output_paths[[i]] <- out$res
          }else if(path_type == 'interlevel'){
            if(length(interlevel_from) == 0){
              output_paths[[i]] <- NULL
              next
            }
            if(graph.type == 'heterogeneous'){
              seq_ids <- igraph::V(het_g)$name[interlevel_from]
              edgelist <- igraph::as_data_frame(het_g)

              to_cell_nodes <- unique(edgelist$to[edgelist$from == seq_ids & edgelist$type == 'interlevel'])
              to_cell_nodes <- which(igraph::V(het_g)$name %in% to_cell_nodes)

              out <- igraph::all_shortest_paths(het_g, from = interlevel_from, to = to_cell_nodes, weights = NA)
              output_paths[[i]] <- out$res
            }
          }

       }
     }
     names(output_paths) <- paths

     return(output_paths)
  }


  prepare_for_plot <- function(tree){

    g <- tree@tree
    paths <- tree@paths
    path_types <- names(paths)[names(paths)!='interlevel']
    df <- igraph::as_data_frame(g, what = 'vertices')

    out_df <- vector(mode = 'list', length = length(path_types))

    for(i in 1:length(path_types)){

      node_ids <- unname(unlist(paths[path_types[i]]))


      if(is.null(color.by)){
        color.by <- 'node_type'
      }

      if(cell.frequency){
        if(color.by == 'node_type'){
          features <- df[color.by][node_ids,]
          counts <- df$cell_number[node_ids]
        }else{
          features <- unlist(df[color.by][node_ids,])
          counts <- unlist(df[paste0(color.by, '_counts')][node_ids,])
        }

      }else{
        if(color.by == 'node_type'){
          features <- df[color.by][node_ids,]
          counts <- rep(1, length(features))
        }else{
          features <- df[color.by][node_ids,]
          counts <- df[paste0(color.by, '_counts')][node_ids,]
          features <- mapply(function(x,y) x[which.max(y)], features, counts)
          counts <- rep(1, length(features))
        }
      }


      out_df[[i]] <- data.frame(features = unlist(features), counts = unlist(counts))
      out_df[[i]] <- out_df[[i]] %>%
                     dplyr::group_by(features) %>%
                     dplyr::summarise(counts = sum(counts))

      out_df[[i]]$sample_id <- paste0(tree@sample_id, collapse = ';')
      out_df[[i]]$clonotype_id <- paste0(tree@clonotype_id, collapse = ';')
      out_df[[i]]$path_type <- path_types[i]
      out_df[[i]]$total_nodes <- length(node_ids)
      out_df[[i]]$total_cells <- sum(df$cell_number[node_ids])
    }
    out_df <- do.call('rbind', out_df)

    return(out_df)
  }


  if(inherits(trees, 'list')){

    for(i in 1:length(trees)){
      final_paths <- vector(mode = "list", length = length(trees[[i]]))

      if(parallel){
        requireNamespace('parallel')
        cores <- parallel::detectCores()
        final_paths <- parallel::mclapply(trees[[i]], mc.cores = cores,
                                                 FUN = function(x) {x %>% get_paths()

                                                                    })

      }else{
        final_paths <- lapply(trees[[i]], function(x) {x %>% get_paths()
                                                              })
      }

      for(j in 1:length(trees[[i]])){
        trees[[i]][[j]]@paths <- final_paths[[j]]
      }



    }

  }else if(inherits(trees, 'AntibodyForests')){
    trees@paths <- trees %>% get_paths()


  }else{
    stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
  }

  if(plot.results){

    if(!inherits(trees, 'list')){
      trees <- list(trees)
    }

    plot_df <- vector(mode = 'list', length = length(trees))

    for(i in 1:length(trees)){
      if(inherits(trees[[i]], 'list')){
        plot_df[[i]] <- lapply(trees[[i]], function(x) prepare_for_plot(x))
      }else{
        plot_df[[i]] <- prepare_for_plot(trees[[i]])
      }
      plot_df[[i]] <- do.call('rbind', plot_df[[i]])

      #Could not find a better solution via dplyr........
      path_types <- unique(plot_df[[i]]$path_type)

      if(cell.frequency){
        plot <-  ggplot2::ggplot(plot_df[[i]], ggplot2::aes(fill = features, y = counts, x = stats::reorder(clonotype_id, -total_cells)))
        plot <- plot + ggplot2::labs(title = unique(plot_df[[i]]$sample_id), x = "Path per clonotype", y = "Total cell counts per path")
      }else{
        plot <-  ggplot2::ggplot(plot_df[[i]], ggplot2::aes(fill = features, y = counts, x = stats::reorder(clonotype_id, -total_nodes)))
        plot <-  plot + ggplot2::labs(title = unique(plot_df[[i]]$sample_id), x = "Path per clonotype", y = "Total node counts per path")
      }

      plot <- plot +
             ggplot2::geom_bar(stat="identity", width=0.6, color="black") +
             ggplot2::scale_y_continuous(expand = c(0,0)) +
             ggplot2::theme_bw() +
             ggplot2::theme_classic() +
             ggplot2::ggtitle(unique(plot_df[[i]]$sample_id)) +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x =  ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))


      if(length(path_types) > 1){
        plot <- plot + ggplot2::facet_wrap(~path_type, scales = "free_x")
      }

      plot(plot)

    }
    plot_df <- do.call('rbind', plot_df)

    if(length(unique(plot_df$sample_id)) > 1){
      if(cell.frequency){
        violin_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample_id, fill = sample_id, y = total_cells)) +
                       ggplot2::labs(x = "Sample id", y = "Cell counts per path")

      }else{
        violin_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample_id, fill = sample_id, y = total_nodes)) +
                       ggplot2::labs(x = "Sample id", y = "Node counts per path")

      }
      violin_plot <- violin_plot +
                     ggplot2::geom_violin(trim = F) +
                     ggplot2::geom_boxplot(width=0.1, fill = 'white') +
                     ggplot2::theme_bw() +
                     ggplot2::theme_classic() +
                     ggplot2::theme(legend.position = "none", text = ggplot2::element_text(size = 20)) +
                     ggplot2::scale_fill_brewer(palette="Dark2")

      plot(violin_plot)
    }

  }

  return(trees)
}
