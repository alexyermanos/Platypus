#' Calculates the node transitions frequencies for a given feature and an AntibodyForests object

#'@description Node transitions represent the number of (un)directed edges between two feature values of a given feature type in a single sequence similarity network or minimum spanning trees (e.g., between VDJ_cgene = 'IGHA' and VDJ_cgene = 'IGHM'). The resulting AntibodyForests objects will contain a new slot - node_transitions. Will also output bar plots of the transition frequencies.

#' @param trees nested list of AntibodyForests objects or single object, as obtained from the AntibodyForests function.
#' @param features vector of strings - features for the node transition counting. Each node will be assigned the dominant feature if its constituent cells have different feature values (e.g., for features = c('seurat_clusters')). These features need to be first added as vertex attributes when creating the AntibodyForests networks, specified in the node.features parameter of AntibodyForests.
#' @param combined boolean - if T, will assign a new feature by combining the feature values as specified in the features parameter (e.g., if node 1 has VDJ_cgene = 'IGHM' and seurat_clusters = 1, will create a new feature = 'IGHM 1' for that node, for joint node transitions).
#' @param graph.type string - the graph type available in the AntibodyForests object which will be used as the function input.
#' Currently supported network/analysis types: 'tree' (for the minimum spanning trees or sequence similarity networks obtained from the main AntibodyForests function), 'heterogeneous' for the bipartite graphs obtained via AntibodyForests_heterogeneous, 'dynamic' for the dynamic networks obtained from AntibodyForests_dynamics.
#' @param permutation.test boolean - if T, will perform a permutation statistical test on the node feature transition values.
#' @param exclude.germline boolean - if T, will exclude the germline from the node feature transitions counting.
#' @param exclude.intermediates boolean - if T, will exclude the intermediate nodes expanded using AntibodyForests_expand_intermediates() from the node transitions counting.
#' @param n.permutations integer - number of node feature permutations for the permutation test.
#' @param plot.results boolean - if T, will display the node transitions counts as bar plots, per feature specified in the features parameter.
#' @param parallel boolean - whether to execute the main subroutine in parallel or not. Requires the 'parallel' R package to be installed.

#' @return nested list of AntibodyForests objects for each clonotype and each sample or single object, with an additional node_transitions slot of transition/edge counts.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_node_transitions(trees,
#' graph.type = 'tree', features = 'VDJ_cgene',
#' plot.results = T)
#'}



AntibodyForests_node_transitions <- function(trees,
                                             features,
                                             combined,
                                             graph.type,
                                             permutation.test,
                                             exclude.germline,
                                             exclude.intermediates,
                                             n.permutations,
                                             plot.results,
                                             parallel){

  if(missing(trees)) stop('Please input a nested list of AntibodyForests objects - the output of AntibodyForests')
  if(missing(features)) features <- NULL
  if(missing(combined)) combined <- F
  if(missing(graph.type)) graph.type <- 'tree'
  if(missing(permutation.test)) permutation.test <- T
  if(missing(exclude.germline)) exclude.germline <- T
  if(missing(exclude.intermediates)) exclude.intermediates <- T
  if(missing(n.permutations)) n.permutations <- 10
  if(missing(plot.results)) plot.results <- T
  if(missing(parallel)) parallel <- T

  feature_name <- NULL
  counts <- NULL
  parent_child <- NULL
  label_both <- NULL

  get_node_transitions <- function(tree){

    if(is.null(features)){
      features <- tree@feature_names

      if(is.null(features)){
        features <- 'node_type'
      }
    }else{
      features <- features
    }

    if(graph.type == 'heterogeneous'){
      g <- tree@heterogeneous
    }else{
      g <- tree@tree
    }

    if(combined){
      feature_values <- ''

      for(feature in features){
        selected_values <- igraph::vertex_attr(g, name=feature)
        counts <- igraph::vertex_attr(g, name=paste0(feature, '_counts'))
        max_indices <- lapply(counts, function(x) which.max(x))
        max_values <- unlist(mapply(function(x,y) x[y], selected_values, max_indices))

        feature_values <- paste0(feature_values, '  ', max_values)
      }

      feature_values <- lapply(feature_values, function(x) substring(x, 3))
      features <- paste0(features, collapse=' ')

      g <- igraph::set_vertex_attr(g, name=features, value=unlist(feature_values))
    }



    edgelist <- igraph::as_edgelist(g)

    if(exclude.intermediates & ('intermediate' %in% igraph::V(g)$node_type)){
      last_node_index <- length(which(igraph::V(g)$node_type=='sequence'))
      if(any(igraph::V(g)$node_type=='germline')){
        last_node_index <- which(igraph::V(g)$node_type=='germline')
      }
      node_to_intermediate <- edgelist[which(edgelist[,1] <= last_node_index & edgelist[,2] > last_node_index),]
      intermediate_to_node <- edgelist[which(edgelist[,1] > last_node_index & edgelist[,2] <= last_node_index),]
      if(is.null(nrow(node_to_intermediate))){
        node_to_intermediate <- t(as.matrix(node_to_intermediate))
      }
      if(is.null(nrow(intermediate_to_node))){
        intermediate_to_node <- t(as.matrix(intermediate_to_node))
      }

      node_to_node <- cbind(node_to_intermediate[,1], intermediate_to_node[,2])
      edgelist <- edgelist[which(edgelist[,1] <= last_node_index & edgelist[,2] <= last_node_index),]
      edgelist <- rbind(edgelist, node_to_node)
      edgelist <- edgelist[order(edgelist[,1], edgelist[,2]),]

    }
    unique_pairs_dfs <- list()
    for(i in 1:length(features)){
      #We need first to exclude list values per node - in the future, might add separate vertex attr called 'OVA_binder_pie' and 'OVA_binder_pie_counts'

      if(length(unlist(igraph::vertex_attr(g, name=features[i]))) > length(igraph::V(g)) ){
        selected_values <- igraph::vertex_attr(g, name=features[i])
        counts <- igraph::vertex_attr(g, name=paste0(features[i], '_counts'))
        max_indices <- lapply(counts, function(x) which.max(x))
        max_values <- unlist(mapply(function(x,y) x[y], selected_values, max_indices))
        g<-igraph::set_vertex_attr(g, name=features[i], value=max_values)
      }


      from_values <- unlist(igraph::vertex_attr(g, name=features[i], index=edgelist[,1]))
      to_values <- unlist(igraph::vertex_attr(g, name=features[i], index=edgelist[,2]))

      df <- data.frame(from=from_values, to=to_values)
      unique_pairs_dfs[[i]] <- df[!duplicated(df),]
      unique_pairs_dfs[[i]] <- unique_pairs_dfs[[i]][order(unique_pairs_dfs[[i]][,1], unique_pairs_dfs[[i]][,2]),]
      row.names(unique_pairs_dfs[[i]]) <- NULL

      pair_counts <- lapply(paste0(unique_pairs_dfs[[i]][,1], ' ', unique_pairs_dfs[[i]][,2]), function(x) length(which(paste0(df[,1], ' ', df[,2])==x)))
      unique_pairs_dfs[[i]]$counts <- unlist(pair_counts)
      unique_pairs_dfs[[i]]$feature_name <- features[i]
    }

    output_df <- do.call('rbind', unique_pairs_dfs)

    ordered_clonotypes <- unique(unlist(igraph::V(g)$clonotype_id[which(igraph::V(g)$node_type=='sequence')]))
    ordered_clonotypes <- ordered_clonotypes[order(nchar(ordered_clonotypes), ordered_clonotypes)]
    ordered_clonotypes <- paste0(ordered_clonotypes, collapse=';')
    output_df$clonotype_id <- rep(ordered_clonotypes, nrow(output_df))


    ordered_samples <- unique(unlist(igraph::V(g)$sample_id[which(igraph::V(g)$node_type=='sequence')]))
    ordered_samples <- ordered_samples[order(nchar(ordered_samples), ordered_samples)]
    ordered_samples <- paste0(ordered_samples, collapse=';')
    output_df$sample_id <- rep(ordered_samples, nrow(output_df))

    if(exclude.germline){
      output_df <-  output_df[which(!(stringr::str_detect(output_df$from, 'germline'))),]
      output_df <-  output_df[which(!(stringr::str_detect(output_df$to, 'germline'))),]
    }

    return(output_df)
  }


  permute_graph_nodes <- function(tree){
    if(is.null(features)){
      features <- tree@feature_names

      if(is.null(features)){
        features <- 'node_type'
      }
    }else{
      features <- features
    }

    if(graph.type == 'heterogeneous'){
      g <- tree@heterogeneous
    }else{
      g <- tree@tree
    }

    if(exclude.germline & exclude.intermediates){
      ids <- which(igraph::V(g)$node_type != 'germline' & igraph::V(g)$node_type != 'intermediate')
    }else if(exclude.germline){
      ids <- which(igraph::V(g)$node_type != 'germline')
    }else if(exclude.intermediates){
      ids <- which(igraph::V(g)$node_type != 'intermediate')
    }else{
      ids <- 1:length(igraph::V(g))
    }


    for(feature in features){
      selected_values <- igraph::vertex_attr(g, name = feature)[ids]
      counts <- igraph::vertex_attr(g, name=paste0(feature, '_counts'))[ids]
      max_indices <- lapply(counts, function(x) which.max(x))
      max_values <- unlist(mapply(function(x,y) x[y], selected_values, max_indices))
      max_values[is.na(max_values)] <- 'unknown'

      g <- igraph::set_vertex_attr(g, name = feature, index = ids, value = sample(max_values))
    }

    if(graph.type == 'heterogeneous'){
      tree@heterogeneous <- g
    }else{
      tree@tree <- g
    }

    return(tree)
  }


  permutation_test <- function(output_df, permut_df){
    features <- unique(output_df$feature_name)

    out_dfs <- vector(mode = 'list', length = length(features))

    for(j in 1:length(features)){
      subset_df <- permut_df[which(permut_df$feature_name == features[j]),]
      out_dfs[[j]] <- output_df[which(output_df$feature_name == features[j]),]

      for(i in 1:nrow(out_dfs[[j]])){
        from <- out_dfs[[j]]$from[i]
        to <- out_dfs[[j]]$to[i]

        counts <- subset_df$counts[subset_df$from == from & subset_df$to == to]

        if(length(counts)== 0) {counts <- 0}
        out_dfs[[j]]$expected_counts[i] <- mean(counts)
        f <- stats::ecdf(counts)
        out_dfs[[j]]$p_value[i] <- 1 - f(out_dfs[[j]]$counts[i])
      }
    }

    output_df <- do.call('rbind', out_dfs)
    return(output_df)
  }



  plot_transitions <- function(tree){

    node_transitions <- tree@node_transitions
    node_transitions$parent_child <- paste0(node_transitions$from, ' -- ', node_transitions$to)

    features <- unique(node_transitions$feature_name)

    sample_id <- paste0(tree@sample_id, collapse = '; ')
    clonotype_id <- tree@clonotype_id

    if(length(clonotype_id) > 5){
      clonotype_id <- 'joined clonotypes'
    }else{
      clonotype_id <- paste0(clonotype_id, collapse = '; ')

    }

    barplot <- ggplot2::ggplot(node_transitions, ggplot2::aes(fill = feature_name, y = counts, x = parent_child)) +
               ggplot2::geom_bar(stat="identity", width=0.6, color="black") +
               ggplot2::scale_y_continuous(expand = c(0,0)) +
               ggplot2::theme_bw() +
               ggplot2::theme_classic() +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x =  ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
               ggplot2::labs(title = paste0(sample_id, ' - ', clonotype_id), x = 'Parent-child transitions', y = 'Number of transitions')

    if(length(features) > 1){
      barplot <- barplot + ggplot2::facet_wrap(~feature_name, scales = "free")
    }

    #TO DO: add directed networks w edge widths scaled by transition counts
    permutation_histogram <- list()
    if(!is.null(tree@permuted_transitions)){
      permuted_transitions <- tree@permuted_transitions
      permuted_transitions$parent_child <- paste0(permuted_transitions$from, ' -- ', permuted_transitions$to)

      for(i in 1:length(features)){

        node_transitions_subset <- node_transitions[node_transitions$feature_name == features[i],]
        permuted_transitions_subset <- permuted_transitions[permuted_transitions$feature_name == features[i],]


        permutation_histogram[[i]] <- ggplot2::ggplot(permuted_transitions_subset, ggplot2::aes(fill = feature_name, x = counts)) +
                                     ggplot2::geom_histogram(bins = 30) +
                                     ggplot2::geom_vline(data = node_transitions_subset, ggplot2::aes(xintercept = counts),
                                                 color="black", linetype=3, size=0.75) +
                                     ggplot2::theme_bw() +
                                     #ggplot2::theme_classic() +
                                     ggplot2::facet_grid(to ~ from, labeller = label_both, scales="free") +
                                     ggplot2::labs(title = paste0(sample_id, ' - ', clonotype_id))

      }
    }


    return(list(barplot = barplot,
                permutation_histogram = permutation_histogram))

  }


  if(inherits(trees, 'list')){

    for(i in 1:length(trees)){
      transitions_dfs <- vector(mode = "list", length = length(trees[[i]]))

      if(parallel){
        requireNamespace('parallel')
        cores <- parallel::detectCores()
        transitions_dfs <- parallel::mclapply(trees[[i]], mc.cores = cores,
                                                 FUN = function(x) {x %>% get_node_transitions()
                                                                    })

      }else{
        transitions_dfs <- lapply(trees[[i]], function(x) {x %>% get_node_transitions()
                                                              })
      }

      for(j in 1:length(trees[[i]])){

        if(permutation.test){
          requireNamespace('parallel')
          cores <- parallel::detectCores()
          permuted_trees <- parallel::mclapply(1:n.permutations, mc.cores = cores, FUN = function(x) {trees[[i]][[j]] %>% permute_graph_nodes() })
          permuted_transitions <-  parallel::mclapply(permuted_trees, mc.cores = cores, FUN = function(x) {x %>% get_node_transitions() })
          permuted_transitions <- do.call('rbind', permuted_transitions)
          trees[[i]][[j]]@permuted_transitions <- permuted_transitions
          transitions_dfs[[j]] <- permutation_test(transitions_dfs[[j]], permuted_transitions)
        }

        trees[[i]][[j]]@node_transitions <- transitions_dfs[[j]]

      }
    }

    if(plot.results){
      for(i in 1:length(trees)){
        for(j in 1:length(trees[[i]])){

          output_plots <- trees[[i]][[j]] %>% plot_transitions()
          plot(output_plots$barplot)
          if(length(output_plots$permutation_histogram) > 0){
            histograms <- plots$permutation_histogram
            for(i in 1:length(histograms)){
              plot(histograms[[i]])

            }
          }
        }
      }
    }

  }else if(inherits(trees, 'AntibodyForests')){
    node_transitions <- trees %>% get_node_transitions()

    if(permutation.test){
      requireNamespace('parallel')
      cores <- parallel::detectCores()
      permuted_trees <- parallel::mclapply(1:n.permutations, mc.cores = cores, FUN = function(x) {trees %>% permute_graph_nodes() })
      permuted_transitions <-  parallel::mclapply(permuted_trees, mc.cores = cores, FUN = function(x) {x %>% get_node_transitions()})
      permuted_transitions <- do.call('rbind', permuted_transitions)
      trees@permuted_transitions <- permuted_transitions
      node_transitions <- permutation_test(node_transitions, permuted_transitions)
    }

    trees@node_transitions <- node_transitions

    if(plot.results){
      plots <- trees %>% plot_transitions()
      plot(plots$barplot)
      if(length(plots$permutation_histogram) > 0){
        histograms <- plots$permutation_histogram
        for(i in 1:length(histograms)){
          plot(histograms[[i]])
        }
      }
    }

  }else{
    stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
  }


  return(trees)

}
