AntibodyForests_communities <- function(trees,
                                        VGM,
                                        community.algorithm,
                                        graph.type,
                                        which.bipartite,
                                        features,
                                        count.level,
                                        additional.parameters){

  if(missing(trees)) stop('Please input a nested list of AntibodyForests objects/ an AntibodyForests objects for community analysis')
  if(missing(VGM)) VGM <- NULL
  if(missing(community.algorithm)) community.algorithm <- 'louvain'
  if(missing(graph.type)) graph.type <- 'tree'
  if(missing(which.bipartite)) which.bipartite <- 'both'
  if(missing(features)) features <- NULL
  if(missing(count.level)) count.level <- 'cells'
  if(missing(additional.parameters)) additional.parameters <- list()

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
        g <- sequence_g
      }else if(which.bipartite == 'cells'){
        g <- cell_g
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

  assemble_bipartite <- function(graphs, original_graph){

    sequence_g <- graphs[[1]]
    cell_g <- graphs[[2]]

    sequence_vertices <- igraph::as_data_frame(sequence_g, what = 'vertices')
    cell_vertices <- igraph::as_data_frame(cell_g, what = 'vertices')

    edgelist <- igraph::as_edgelist(original_graph)
    g <- igraph::graph_from_data_frame(d = edgelist, directed = F, vertices = rbind(sequence_vertices, cell_vertices))

    return(g)

  }

  community_detection <- function(g, original_graph){
    bipartite_type <- names(g)
    communities <- vector(mode = 'list', length = length(g))

    for(i in 1:length(g)){
      if(bipartite_type[i] == 'sequence'){
        weights <- igraph::E(g[[i]])$weight
        igraph::E(g[[i]])$weight <- 1/weights
      }

      if(community.algorithm[i] == 'edge_betweenness'){
        community <- rlang::exec(igraph::cluster_edge_betweenness, g[[i]], !!!additional.parameters)$membership

      }else if(community.algorithm[i] == 'fast_greedy'){
        community <- rlang::exec(igraph::cluster_fast_greedy, g[[i]], !!!additional.parameters)$membership

      }else if(community.algorithm[i] == 'label_prop'){
        community <- rlang::exec(igraph::cluster_label_prop, g[[i]], !!!additional.parameters)$membership

      }else if(community.algorithm[i] == 'leading_eigen'){
        community <- rlang::exec(igraph::cluster_leading_eigen, g[[i]], !!!additional.parameters)$membership

      }else if(community.algorithm[i] == 'louvain'){
        community <- rlang::exec(igraph::cluster_louvain, g[[i]], !!!additional.parameters)$membership

      }else if(community.algorithm[i] == 'leiden'){
        community <- rlang::exec(igraph::cluster_leiden, g[[i]], !!!additional.parameters)$membership

      }else if(community.algorithm[i] == 'optimal'){
        community <- rlang::exec(igraph::cluster_optimal, g[[i]], !!!additional.parameters)$membership

      }else if(community.algorithm[i] == 'spinglass'){
        community <- rlang::exec(igraph::cluster_spinglass, g[[i]], !!!additional.parameters)$membership

      }else if(community.algorithm[i] == 'walktrap'){
        community <- rlang::exec(igraph::cluster_walktrap, g[[i]], !!!additional.parameters)$membership

      }else{
        stop('Unrecognized community detection algorithm!')
      }

      if(bipartite_type[i] == 'cell'){
        community <- paste0(community, '_', bipartite_type[i])
      }

      g[[i]] <- igraph::set_vertex_attr(g[[i]], name = 'community', value = as.character(community))

      if(bipartite_type[i] == 'sequence'){
        igraph::E(g[[i]])$weight <- weights
      }

    }

    if(length(g) == 2){
      g <- assemble_bipartite(g, original_graph)
    }else{
      g <- g[[1]]
    }

    return(g)
  }


  community_barplot <- function(g){

    vertex_df <- igraph::as_data_frame(g, what = 'vertices')
    clusters <- unique(vertex_df$community)

    final_dfs <- vector(mode = 'list', length = length(features))
    for(i in 1:length(features)){

      feat <- features[i]
      cluster_dfs <- vector(mode = 'list', length = length(clusters))

      for(j in 1:length(clusters)){
        cluster <- clusters[j]
        df <- vertex_df[vertex_df$community == cluster,]

        if(count.level == 'cells'){
          if(paste0(feat, '_counts') %in% colnames(df)){
            feature_values <- unlist(df[[feat]])
            counts <- unlist(df[[paste0(feat, '_counts')]])

          }else{
            feature_values <- unlist(df[[feat]])
            counts <- df[['cell_number']]
          }

        }else{
          if(paste0(feat, '_counts') %in% colnames(df)){
            counts <- df[[paste0(feat, '_counts')]]
            feature_values <- df[[feat]]
            feature_values <- mapply(function(x, y) x[which.max(y)], feature_values, counts)
            counts <- rep(1, length(feature_values))

          }else{
            feature_values <- unlist(df[[feat]])
            counts <- rep(1, length(feature_values))

          }
        }

        feature_values[is.na(feature_values) | feature_values == ''] <- 'unknown'
        unique_features <- unique(feature_values)
        counts_per_feature <- sapply(unique_features, function(x) sum(counts[which(feature_values == x)]))

        cluster_dfs[[j]] <- data.frame(features = unlist(unique_features), counts = unlist(unname(counts_per_feature)), feature_name = feat, cluster = cluster)
      }

      final_dfs[[i]] <- do.call('rbind', cluster_dfs)

    }

    final_dfs <- do.call('rbind', final_dfs)

    out_plots <- list()
    for(i in 1:length(features)){
      df_subset <- final_dfs[final_dfs$feature_name == features[i],]

      out_plots <- ggplot2::ggplot(df_subset, ggplot2::aes(fill = features, y = counts, x = reorder(cluster, nchar(cluster)))) +
                   ggplot2::geom_bar(stat="identity", width=0.6, color="black") +
                   ggplot2::scale_y_continuous(expand = c(0,0)) +
                   ggplot2::theme_bw() +
                   ggplot2::theme_classic() +
                   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x =  ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
                   ggplot2::labs(title = features[i])

      if(count.level == 'cells'){
        out_plots <- out_plots + ggplot2::labs(x = "Cluster", y = "Number of cells")
      }else{
        out_plots <- out_plots + ggplot2::labs(x = "Cluster", y = "Number of nodes")
      }

    }

    return(out_plots)
  }

  append_community <- function(g, VGM){
    barcodes <- igraph::V(g)$cell_barcodes
    clusters <- igraph::V(g)$community

    clusters <- mapply(function(x,y) rep(y, length(x)), barcodes, clusters)

    barcode_df <- data.frame(cbind(unlist(barcodes), unlist(clusters)))
    colnames(barcode_df) <- c('barcode', 'community')

    if(!is.null(VGM[[1]])){
      VGM[[1]] <- merge(VGM[[1]], barcode_df, by = 'barcode', all.x = T)
    }

    if(!is.null(VGM[[2]])){
      metadata <- VGM[[2]]@meta.data
      metadata$barcode <- rownames(metadata)
      metadata <- merge(metadata, barcode_df, by = 'barcode', all.x = T)
      VGM[[2]] <- Seurat::AddMetaData(object = VGM[[2]], metadata = as.factor(metadata$community), col.name = 'community')
    }

    return(VGM)
  }

  features <- get_feature_names(trees, features)

  if(inherits(trees, 'list')){
    for(i in 1:length(trees)){
      for(j in 1:length(trees[[i]])){
        if(graph.type == 'tree'){
          g <- trees[[i]][[j]] %>% get_graph() %>% community_detection(trees[[i]][[j]]@tree)

          if(!is.null(VGM)){
            output_vgm <- append_community(g, VGM)
          }

          trees[[i]][[j]]@tree <- g
          p <- community_barplot(g)
          plot(p)


        }else if(graph.type == 'dynamic'){
          g <- trees[[i]][[j]] %>% get_graph() %>% community_detection(trees[[i]][[j]]@dynamic)

          if(!is.null(VGM)){
            output_vgm <- append_community(g, VGM)
          }

          trees[[i]][[j]]@dynamic <- g
          p <- community_barplot(g)
          plot(p)


        }else if(graph.type == 'heterogeneous'){
          g <- trees[[i]][[j]] %>% get_graph() %>% community_detection(trees[[i]][[j]]@heterogeneous)

          if(!is.null(VGM)){
            output_vgm <- append_community(g, VGM)
          }

          trees[[i]][[j]]@heterogeneous <- g
          p <- community_barplot(g)
          plot(p)


        }else{
          stop('Unrecognized graph type!')
        }

      }
    }

  }else if(inherits(trees, 'AntibodyForests')){
    if(graph.type == 'tree'){
      g <- trees %>% get_graph() %>% community_detection(trees@tree)

      if(!is.null(VGM)){
        output_vgm <- append_community(g, VGM)
      }

      trees@tree <- g
      p <- community_barplot(g)
      plot(p)


    }else if(graph.type == 'dynamic'){
      g <- trees %>% get_graph() %>% community_detection(trees@dynamic)

      if(!is.null(VGM)){
        output_vgm <- append_community(g, VGM)
      }

      trees@dynamic <- g
      p <- community_barplot(g)
      plot(p)


    }else if(graph.type == 'heterogeneous'){
      g <- trees %>% get_graph() %>% community_detection(trees@heterogeneous)

      if(!is.null(VGM)){
        output_vgm <- append_community(g, VGM)
      }

      trees@heterogeneous <- g
      p <- community_barplot(g)
      plot(p)


    }else{
      stop('Unrecognized graph type!')
    }


  }else{
    stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
  }


  if(!is.null(VGM)){
    return(output_vgm)
  }else{
    return(trees)
  }
}
