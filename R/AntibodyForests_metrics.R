
#Further analysis for bipartite - metacells on the cell graph, overlaps/all that in the resulting incidence matrix (now we ensure some overlap btw sequences and their corresponding cells)
#add incidence matrix analyses for metacell bipartite graphs

AntibodyForests_metrics <- function(trees,
                                    metrics,
                                    graph.type,
                                    features,
                                    exclude.intermediates,
                                    exclude.germline,
                                    separate.bipartite,
                                    parallel
                                    ){


  if(missing(trees)) stop('Please input your nested network list, obtained from AntibodyForests_parallel')
  if(missing(metrics)) metrics <- c('node_metrics')
  if(missing(graph.type)) graph.type <- 'tree'
  if(missing(features)) features <- NULL
  if(missing(exclude.intermediates)) exclude.intermediates <- T
  if(missing(exclude.germline)) exclude.germline <- T
  if(missing(separate.bipartite)) separate.bipartite <- F
  if(missing(parallel)) parallel <- T


  node_metrics <- c('betweenness', 'closeness', 'eigenvector', 'authority_score',
                    'local_cluster_coefficient', 'average_cluster_coefficient',
                    'strength', 'degree',
                    'eccentricity', 'pagerank', 'daughters',
                    'path_from_germline', 'weighted_path_from_germline',
                    'path_from_most_expanded',
                    'weighted_path_from_most_expanded',
                    'path_from_hub', 'weighted_path_from_hub')


  graph_metrics <- c('average_expansion', 'average_daughters', 'average_degree', 'degree_centrality', 'linear',
                     'dominant_features', 'feature_counts_nodes')


  if('all' %in% metrics){
    metrics <- c(node_metrics, graph_metrics)
  }else if('node_metrics' %in% metrics){
    metrics <- node_metrics
  }else if('graph_metrics' %in% metrics){
    metrics <- graph_metrics
   }


  get_graph <- function(tree){

    if(graph.type == 'tree'){
      g <- tree@tree

    }else if(graph.type == 'heterogeneous'){
      g <- tree@heterogeneous

      if(separate.bipartite){
        cell_vertices <- which(igraph::V(g)$type == 'cell')
        sequence_vertices <- which(igraph::V(g)$type == 'sequence')
        sequence_g <- igraph::delete_vertices(g, cell_vertices)
        cell_g <- igraph::delete_vertices(g, sequence_vertices)

        g <- list()
        g[[1]] <- sequence_g
        g[[2]] <- cell_g
        names(g) <- c('sequence', 'cell')
      }

    }else if(graph.type == 'dynamic'){
      g <- tree@dynamic

    }else{
      stop('Graph type not found!')
    }

    if(is.null(g)){
      stop(paste0('Could not find the ', graph.type, ' graph!'))
    }

    return(g)
  }


  compute_node_metrics <- function(g, feature_names){

    if('betweenness' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'betweenness' , value = igraph::betweenness(g, weights = NULL))
    }

    if('closeness' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'closeness' , value = igraph::closeness(g, mode = c("all"), weights = NULL, normalized = TRUE))
    }

    if('eigenvector' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'eigenvector' , value = igraph::eigen_centrality(g, scale=T, weights=NULL)$vector)
    }

    if('authority_score' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'authority_score' , value = unname(igraph::authority_score(g)$vector))
    }

    if('degree' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'degree' , value = igraph::degree(g))
    }

    if('strength' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'strength' , value = igraph::strength(g))
    }

    if('local_cluster_coefficient' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'local_cluster_coefficient', value = igraph::transitivity(g, type = "local", vids=NULL))
    }

    if('average_cluster_coefficient' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'average_cluster_coefficient', value = igraph::transitivity(g, type = "average"))
    }

    if('daughters' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'daughters', value = igraph::degree(g, mode='out'))
    }

    if('pagerank' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'pagerank' , value = igraph::page_rank(g, directed=T, weights=NULL)$vector)
    }

    if('eccentricity' %in% metrics){
      g <- igraph::set_vertex_attr(g, name = 'eccentricity', value = igraph::eccentricity(g))
    }


    if('path_from_germline' %in% metrics | 'weighted_path_from_germline' %in% metrics){

      root_index = which(igraph::V(g)$node_type=='germline')

      #TO DO : GENERALIZE THE DISTANCE CALCULATION BY PICKING ROOTS BEFORE (OR KEEP LIKE THIS BCS MORE VISIBLE/EASIER TO UNDERSTAND)
      if('path_from_germline' %in% metrics){
        if(length(root_index) == 0){
          igraph::V(g)$path_from_germline <- rep(NA, length(igraph::V(g)))

        }else if(length(root_index) > 1 ){
          #joined trees case
          for(root in rev(root_index)){
            distances <- igraph::distances(g, v=root, weights=NA)
            non_inf_distances <- distances[distances!='Inf']
            igraph::V(g)$path_from_germline[(root-length(non_inf_distances)+1):root] <- non_inf_distances
          }
        }else{
          distances <- igraph::distances(g, v=root_index, weights=NA)
          igraph::V(g)$path_from_germline <- distances
        }
      }

      if('weighted_path_from_germline' %in% metrics){
       if(length(root_index) == 0){
         igraph::V(g)$weighted_path_from_germline <- rep(NA, length(igraph::V(g)))

       }else if(length(root_index) > 1 ){ #joined trees case
          for(root in rev(root_index)){
            distances <- igraph::distances(g, v=root, weights=NULL)
            non_inf_distances <- distances[distances!='Inf']
            igraph::V(g)$weighted_path_from_germline[(root-length(non_inf_distances)+1):root] <- non_inf_distances
          }
        }else{
          distances <- igraph::distances(g, v=root_index, weights=NULL)
          igraph::V(g)$weighted_path_from_germline <- distances
        }
      }
    }

      if('path_from_most_expanded' %in% metrics | 'weighted_path_from_most_expanded' %in% metrics){
        root_index <- which(igraph::V(g)$most_expanded=='yes')
        if('path_from_most_expanded' %in% metrics){
          if(length(root_index) == 0){
            igraph::V(g)$path_from_most_expanded <- rep(NA, length(igraph::V(g)))

          }else if(length(root_index) > 1){
            for(root in rev(root_index)){
              distances <- igraph::distances(g, v=root, weights=NA)
              non_inf_distances <- distances[distances!='Inf']
              igraph::V(g)$path_from_most_expanded[root:(root+length(non_inf_distances)-1)] <- non_inf_distances #Priors which might break: most_expandeds are always the first sequence in a unique clonotype, germlines are always the last
            }
          }else{

            distances <- igraph::distances(g, v=root_index, weights=NA)
            igraph::V(g)$path_from_most_expanded <- distances #Priors which might break: most_expandeds are always the first sequence in a unique clonotype, germlines are always the last
          }
        }

        if('weighted_path_from_most_expanded' %in% metrics){

          if(length(root_index) == 0){
            igraph::V(g)$weighted_path_from_most_expanded <- rep(NA, length(igraph::V(g)))

          }else if(length(root_index) > 1){
            for(root in rev(root_index)){
              distances <- igraph::distances(g, v=root, weights=NULL)
              non_inf_distances <- distances[distances!='Inf']
              igraph::V(g)$weighted_path_from_most_expanded[root:(root+length(non_inf_distances)-1)] <- non_inf_distances #Priors which might break: most_expandeds are always the first sequence in a unique clonotype, germlines are always the last
            }
          }else{
            distances <- igraph::distances(g, v=root_index, weights=NULL)
            igraph::V(g)$weighted_path_from_most_expanded <- distances #Priors which might break: most_expandeds are always the first sequence in a unique clonotype, germlines are always the last
          }
        }
      }

      if('path_from_hub' %in% metrics | 'weighted_path_from_hub' %in% metrics){
        root_index <- which(igraph::V(g)$hub=='yes')
        if('path_from_hub' %in% metrics){
          if(length(root_index) == 0){
            igraph::V(g)$path_from_hub <- rep(NA, length(igraph::V(g)))

          }else if(length(root_index) > 1){
            for(root in rev(root_index)){
              distances <- igraph::distances(g, v=root, weights=NA)
              non_inf_distances <- distances[distances!='Inf']
              sample <- igraph::V(g)$sample_id[root]
              clonotype <- igraph::V(g)$clonotype_id[root]
              igraph::V(g)$path_from_hub[igraph::V(g)$sample_id==sample & igraph::V(g)$clonotype_id==clonotype] <- non_inf_distances
            }
          }else{
            distances <- igraph::distances(g, v=root_index, weights=NA)
            igraph::V(g)$path_from_hub <- distances
          }
        }

        if('weighted_path_from_hub' %in% metrics){
          if(length(root_index) == 0){
            igraph::V(g)$weighted_path_from_hub <- rep(NA, length(igraph::V(g)))

          }else if(length(root_index) > 1){
            for(root in rev(root_index)){
              distances <- igraph::distances(g, v=root, weights=NULL)
              non_inf_distances <- distances[distances!='Inf']
              sample <- igraph::V(g)$sample_id[root]
              clonotype <- igraph::V(g)$clonotype_id[root]
              igraph::V(g)$weighted_path_from_hub[igraph::V(g)$sample_id==sample & igraph::V(g)$clonotype_id==clonotype] <- non_inf_distances
            }
          }else{
            distances <- igraph::distances(g, v=root_index, weights=NULL)
            igraph::V(g)$weighted_path_from_hub <- distances
          }
        }
      }

    node_df <- igraph::as_data_frame(g, what='vertices')

    #Postprocessing the node_df dataframe - remove node attributes that are lists (e.g., two diff cell transcriptomic cluster per same sequence) and replace them by the max feature values specified by feature_counts, subsequently remove the feature_counts or replace by cell_numbers
    for(col in feature_names){
      max_indices <- node_df[paste0(col, '_counts')][which(node_df$node_type=='sequence'),]
      max_indices <- unlist(lapply(max_indices, function(x) which.max(x)))
      max_values <- unlist(mapply(function(x,y) x[y], node_df[col][which(node_df$node_type=='sequence'),], max_indices))
      node_df[col][which(node_df$node_type=='sequence'),] <- max_values
      node_df[col] <- unlist(node_df[col]) #to also change the column class from list into numeric/character

      node_df <- node_df[,!(names(node_df) %in% paste0(col, '_counts'))]
    }

    if(exclude.germline){
      node_df <- node_df[node_df$node_type != 'germline',]
    }

    if(exclude.intermediates){
      node_df <- node_df[node_df$node_type != 'intermediate',]

    }


    return(list(graph = g, metrics_df = node_df))
  }


  assemble_bipartite <- function(graphs, original_graph){
    sequence_g <- graphs[[1]]
    cell_g <- graphs[[2]]

    sequence_vertices <- igraph::as_data_frame(sequence_g, what = 'vertices')
    cell_vertices <- igraph::as_data_frame(cell_g, what = 'vertices')

    edgelist <- igraph::as_edgelist(original_graph)
    g <- igraph::graph_from_data_frame(d = edgelist, directed = F, vertices = rbind(sequence_vertices, cell_vertices))
    return(g)
  }


  compute_node_metrics_lapply <- function(tree){

    graphs <- get_graph(tree)
    feature_names <- tree@feature_names

    if(inherits(graphs, 'list')){
      final_graphs <- list()
      metrics_df <- list()

      for(i in 1:length(graphs)){
        out <- compute_node_metrics(g = graphs[[i]], feature_names = feature_names)
        final_graphs[[i]] <- out$graph
        metrics_df[[i]] <- out$metrics_df
      }

      g <- assemble_bipartite(final_graphs, original_graph = tree@heterogeneous)
      metrics_df <- do.call('rbind', metrics_df)


    }else{
      out <- compute_node_metrics(g = graphs, feature_name = feature_names)
      g <- out$graph
      metrics_df <- out$metrics_df
    }

    tree@metrics <- metrics_df

    if(graph.type == 'heterogeneous'){
      tree@heterogeneous <- g
    }else if(graph.type == 'tree'){
      tree@tree <- g
    }else if(graph.type == 'dynamic'){
      tree@dynamic <- g
    }

    return(tree)
  }

  compute_graph_metrics <- function(tree){

    if(!any(metrics %in% graph_metrics)){
      return(tree)
    }

    g <- get_graph(tree)
    #Currently not supported for bipartite/heterogeneous graphs - will add in the future (TO DO: meaningful whole graph metrics for heterogeneous networks)
    if(inherits(g, 'list')){
      g <- g[[1]]
    }

    features <- tree@feature_names

    if(is.null(features)){
      features <- 'node_type'
    }

    sample_id <- paste0(tree@sample_id, collapse = '; ')
    clonotype_id <- tree@clonotype_id

    if(length(clonotype_id) > 5){
      clonotype_id <- 'joined clonotypes'
    }else{
      clonotype_id <- paste0(clonotype_id, collapse = '; ')

    }

    node_df <- igraph::as_data_frame(g, what='vertices')


    #Not subsetting by sequence nodes anymore!
    sequence_subset <- node_df

    if(exclude.intermediates){
      sequence_subset <- node_df[which(node_df$node_type!='intermediate'),]
    }

    if(exclude.germline){
      sequence_subset <- node_df[which(node_df$node_type!='germline'),]
    }

    cell_numbers_df <- NULL
    node_numbers_df <- NULL
    final_graph_df <- NULL

    if(('feature_counts_nodes' %in% metrics) | ('feature_percentages_nodes' %in% metrics)){
      if(is.null(features)){
        stop('Please input a node feature from your igraph objects to calculate counts/percentages')
      }
      total_feature_values <- c()
      total_feature_counts <- c()
      total_feature_names <- c()
      total_feature_percentages <- c()

      for(feature in features){
        feature_values <- unname(unlist(sequence_subset[feature]))
        feature_values[is.na(feature_values)] <- 'NA'
        unique_feature_values <- unlist(unique(feature_values))
        feature_counts <- unlist(lapply(unique_feature_values, function(x) length(which(feature_values==x))))
        feature_percentages <- unlist(feature_counts) / sum(unlist(feature_counts))

        total_feature_values <- c(total_feature_values, unique_feature_values)
        total_feature_counts <- c(total_feature_counts, feature_counts)
        total_feature_names <- c(total_feature_names, rep(feature, length(unique_feature_values)))
        total_feature_percentages <- c(total_feature_percentages, feature_percentages)
      }

      node_numbers_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id,
                                   feature_values=total_feature_values,
                                   feature_counts=total_feature_counts,
                                   feature_percentages=total_feature_percentages,
                                   feature_name=total_feature_names)
    }

    if(('feature_counts_cells' %in% metrics) | ('feature_percentages_cells' %in% metrics)){
      if(is.null(features)){
        stop('Please input a node feature from your igraph objects to calculate counts/percentages')
      }
      total_feature_values <- c()
      total_feature_counts <- c()
      total_feature_names <- c()
      total_feature_percentages <- c()

      for(feature in features){
        feature_values <- unname(unlist(sequence_subset[feature]))
        feature_values[is.na(feature_values)] <- 'NA'
        feature_counts <- unlist(sequence_subset[paste0(feature, '_counts')])

        unique_feature_values <- unlist(unique(feature_values))
        unique_feature_counts <- unlist(lapply(unique_feature_values, function(x) sum(feature_counts[which(feature_values==x)])))

        feature_percentages <- unique_feature_counts / sum(unique_feature_counts)

        total_feature_values <- c(total_feature_values, unique_feature_values)
        total_feature_counts <- c(total_feature_counts, unique_feature_counts)
        total_feature_names <- c(total_feature_names, rep(feature, length(unique_feature_values)))
        total_feature_percentages <- c(total_feature_percentages, feature_percentages)
      }

      cell_numbers_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id,
                                   feature_values=total_feature_values,
                                   feature_counts=total_feature_counts,
                                   feature_percentages=total_feature_percentages,
                                   feature_name=total_feature_names)
    }
    final_graph_df <- rbind(node_numbers_df, cell_numbers_df)


    if('cell_number' %in% metrics){
      cell_number <- sum(unlist(sequence_subset$cell_number[!is.na(sequence_subset$cell_number)]))
      if(is.null(final_graph_df)){
        final_graph_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id, cell_number=cell_number)
      }else{
        final_graph_df$cell_number <- rep(cell_number, nrow(final_graph_df))
      }
    }

    if('average_cell_number' %in% metrics){
      cell_number <- sum(unlist(sequence_subset$cell_number[!is.na(sequence_subset$cell_number)]))
      average_cell_number <- cell_number / nrow(sequence_subset)
      if(is.null(final_graph_df)){
        final_graph_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id, node_number=node_number)
      }else{
        final_graph_df$node_number <- rep(node_number, nrow(final_graph_df))
      }
    }

    if('node_number' %in% metrics){
      node_number <- nrow(sequence_subset)
      if(is.null(final_graph_df)){
        final_graph_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id, node_number=node_number)
      }else{
        final_graph_df$node_number <- rep(node_number, nrow(final_graph_df))
      }
    }

    if('average_expansion' %in% metrics){
      average_expansion <- sum(unlist(sequence_subset$cell_number)) / nrow(sequence_subset)

      if(is.null(final_graph_df)){
        final_graph_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id, average_expansion=average_expansion)
      }else{
        final_graph_df$average_expansion <- rep(average_expansion, nrow(final_graph_df))
      }
    }

    if('average_daughters' %in% metrics){
      indices <- unlist(sequence_subset$label[which(!is.na(sequence_subset$label))])

      daughters <- unlist(igraph::degree(g, mode='out')[indices])
      average_daughters <- sum(daughters) / length(indices)

      if(is.null(final_graph_df)){
        final_graph_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id, average_daughters=average_daughters)
      }else{
        final_graph_df$average_daughters <- rep(average_daughters, nrow(final_graph_df))
      }
    }

    if('average_degree' %in% metrics){
      indices = unlist(sequence_subset$label[which(!is.na(sequence_subset$label))])

      degree <- igraph::degree(g, mode='all')[indices]
      average_degree <- sum(degree) / length(indices)

      if(is.null(final_graph_df)){
        final_graph_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id, average_degree=average_degree)
      }else{
        final_graph_df$average_degree <- rep(average_degree, nrow(final_graph_df))
      }
    }

    if('dominant_features' %in% metrics){

       for(i in 1:length(features)){
         feature_values <- unname(unlist(sequence_subset[feature]))
         feature_values[is.na(feature_values)] <- 'NA'
         unique_feature_values <- unique(feature_values)
         feature_counts <- unlist(lapply(unique_feature_values, function(x) length(which(feature_values==x))))

         max_feature <- unique_feature_values[which.max(feature_counts)]

         if(is.null(final_graph_df)){
           final_graph_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id)
           final_graph_df[paste0('dominant', '_', feature)] <- max_feature
         }else{
           final_graph_df[paste0('dominant', '_', feature)] <- rep(max_feature, nrow(final_graph_df))
         }
       }
      }

    if('linear' %in% metrics){
      indices <- as.vector(sequence_subset$label[which(!is.na(sequence_subset$label))])

      if(!(igraph::is_directed(g))){
        if(any(igraph::degree(g, mode='all')[indices] > 2)){
          linear <- 'no'
        }else{
          linear <- 'yes'
        }
      }else{
        if(any(igraph::degree(g, mode='out')[indices] > 1)){
          linear <- 'no'
        }else{
          linear <- 'yes'
        }
      }

      if(is.null(final_graph_df)){
        final_graph_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id)
        final_graph_df$linear <- linear
      }else{
        final_graph_df$linear <- rep(linear, nrow(final_graph_df))
      }
    }


    if('degree_centrality' %in% metrics){
      degree_centrality <- igraph::centr_degree(g)$centralization

      if(is.null(final_graph_df)){
        final_graph_df <- data.frame(sample_id = sample_id, clonotype_id = clonotype_id)
        final_graph_df$degree_centrality <- degree_centrality
      }else{
        final_graph_df$degree_centrality <- rep(degree_centrality, nrow(final_graph_df))
      }
    }

    if(!is.null(final_graph_df)){
      ordered_clonotypes <- unique(unlist(igraph::V(g)$clonotype_id[which(igraph::V(g)$node_type=='sequence')]))
      ordered_clonotypes <- ordered_clonotypes[order(nchar(ordered_clonotypes), ordered_clonotypes)]
      ordered_clonotypes <- paste0(ordered_clonotypes, collapse=';')
      final_graph_df$clonotype_id <- ordered_clonotypes


      ordered_samples <- unique(unlist(igraph::V(g)$sample_id[which(igraph::V(g)$node_type=='sequence')]))
      ordered_samples <- ordered_samples[order(nchar(ordered_samples), ordered_samples)]
      ordered_samples <- paste0(ordered_samples, collapse=';')
      final_graph_df$sample_id <-ordered_samples
    }

    tree@metrics <- list(tree@metrics, final_graph_df)
    names(tree@metrics) <- c('node_metrics', 'graph_metrics')

    return(tree)
  }

  if(inherits(trees, 'list')){

    for(i in 1:length(trees)){

      if(parallel){
        requireNamespace('parallel')
        cores <- parallel::detectCores()
        trees[[i]] <- parallel::mclapply(trees[[i]], mc.cores = cores, FUN = function(x) x %>% compute_node_metrics_lapply())
        trees[[i]] <- parallel::mclapply(trees[[i]], mc.cores = cores, FUN = function(x) x %>% compute_graph_metrics())



      }else{
        trees[[i]] <- lapply(trees[[i]], function(x)  x %>% compute_node_metrics_lapply())
        trees[[i]] <- lapply(trees[[i]], function(x)  x %>% compute_graph_metrics())


      }
    }

  }else if(inherits(trees, 'AntibodyForests')){
    trees <- trees %>% compute_node_metrics_lapply()
    trees <- trees %>% compute_graph_metrics()

  }else{
    stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
  }


  return(trees)
}
