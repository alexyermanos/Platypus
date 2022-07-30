#' Edge overlap heatmaps for a set of AntibodyForests sequence similarity networks or minimum spanning trees.

#'@description Similar to the AntibodyForests_node_transitions function, will calculate the incidence of features across undirected edges. In this case, each edge will be considered a unique species - with incidence counts across each unique feature value if a specific edge is connected to a node with that feature. Overlap metrics are then calculated for this edge-feature incidence matrix.
#' @param trees nested list of AntibodyForests objects, as obtained from the AntibodyForests function.
#' @param group.by vector of strings - node features to group the edges by (counts edge incidence across the unique feature values for the specified node feature).
#' @param method string - overlap calculator: 'overlap' for unique/public edge counts across the feature values, 'jaccard' to calculate the Jaccard index.

#' @return Edge overlap heatmaps for the specific overlap metric/method.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_overlap(trees, group.by = 'seurat_clusters', method = 'jaccard')
#'}


AntibodyForests_overlap <- function(trees,
                                    group.by,
                                    method
                                    ){
  if(missing(trees)) stop('Please input a nested list of AntibodyForests objects/ a single AntibodyForests object!')
  if(missing(group.by)) group.by <- 'clonotype_id'
  if(missing(method)) method <- 'overlap'

  metric <- NULL

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

  get_edge_incidence <- function(g, features = group.by){
    if(igraph::is_directed(g)){
      g <- igraph::as.undirected(g)
    }

    node_df <- igraph::as_data_frame(g, what = 'vertices')
    edge_df <- igraph::as_data_frame(g, what = 'edges')
    edge_df$edge_id <- 1:nrow(edge_df)

    incidence_dfs <- vector(mode = 'list', length = length(features))

    for(j in 1:length(features)){
      feature <- features[j]

      if(paste0(feature, '_counts') %in% colnames(node_df)){
        feature_values <- node_df[[feature]]
        counts <- node_df[[paste0(feature, '_counts')]]
        max_features <- mapply(function(x,y) x[which.max(y)], feature_values, counts)
        node_df[[feature]] <- max_features
      }

      edge_df$feature_from <- node_df[[feature]][edge_df$from]
      edge_df$feature_from[is.na(edge_df$feature_from) | edge_df$feature_from == ''] <- 'unknown'
      edge_df$feature_to <- node_df[[feature]][edge_df$to]
      edge_df$feature_to[is.na(edge_df$feature_to) | edge_df$feature_to == ''] <- 'unknown'


      unique_features <- unique(c(edge_df$feature_from, edge_df$feature_to))
      unique_features <- unique_features[order(unique_features, nchar(unique_features))]

      incidence_df <- matrix(0, length(unique_features), nrow(edge_df))
      rownames(incidence_df) <- unique_features
      colnames(incidence_df) <- edge_df$edge_id

      for(i in 1:nrow(edge_df)){
        incidence_df[edge_df[i,]$feature_from, edge_df[i,]$edge_id] <- incidence_df[edge_df[i,]$feature_from, edge_df[i,]$edge_id] + 1
        incidence_df[edge_df[i,]$feature_to, edge_df[i,]$edge_id] <- incidence_df[edge_df[i,]$feature_to, edge_df[i,]$edge_id] + 1

      }

      incidence_dfs[[j]] <- incidence_df
    }

    names(incidence_dfs) <- features
    return(incidence_dfs)
  }


  edge_overlap <- function(incidence_dfs){
    final_dfs <- vector(mode = 'list', length(incidence_dfs))

    for(j in 1:length(incidence_dfs)){

      df <- incidence_dfs[[j]]
      combs <- data.frame(t(utils::combn(rownames(df), m = 2, simplify = TRUE)))
      additional <- cbind(rownames(df), rownames(df))
      colnames(additional) <- c('X1', 'X2')
      combs <- rbind(combs, additional)
      combs <- combs[order(combs$X1, combs$X2),]
      combs$metric <- NA
      combs$metric_name <- 'Overlap'
      combs$feature_name <- names(incidence_dfs)[j]


      for(i in 1:nrow(combs)){
        if(combs[i,1] == combs[i,2]){
          combs$metric[i] <- length(which(df[combs[i,1],] == 2))

        }else{
          a <- unname(which(df[combs[i,1],] != 0))
          b <- unname(which(df[combs[i,2],] != 0))

          intersection <- length(intersect(a, b))
          combs$metric[i] <- intersection
        }
      }

      final_dfs[[j]] <- combs
    }

    final_df <- do.call('rbind', final_dfs)
    return(final_df)
  }


  edge_jaccard_similarity <- function(incidence_dfs){

    final_dfs <- vector(mode = 'list', length(incidence_dfs))

    for(j in 1:length(incidence_dfs)){

      df <- incidence_dfs[[j]]
      combs <- data.frame(t(utils::combn(rownames(df), m = 2, simplify = TRUE)))
      combs$metric <- NA
      combs$metric_name <- 'Jaccard similarity'
      combs$feature_name <- names(incidence_dfs)[j]

      for(i in 1:nrow(combs)){
        a <- unname(which(df[combs[i,1],] != 0))
        b <- unname(which(df[combs[i,2],] != 0))

        intersection <- length(intersect(a, b))
        union <- length(a) + length(b) - intersection
        combs$metric[i] <- intersection / union
      }

      final_dfs[[j]] <- combs

    }

    final_df <- do.call('rbind', final_dfs)

    return(final_df)
  }

  plot_overlap <- function(metric_df){

    title_out <- unique(metric_df$metric_name)
    metric_df$metric <- as.numeric(round(metric_df$metric, 3))

    plot_out <- ggplot2::ggplot(metric_df, ggplot2::aes(x = metric_df[,2], y = metric_df[,1], fill = as.numeric(metric))) +
                ggplot2::geom_tile() +
                ggplot2::geom_text(ggplot2::aes(label = metric), size = 4)+
                ggplot2::theme(panel.background = ggplot2::element_blank(),
                            axis.text = ggplot2::element_text(size = 30), axis.line.x = ggplot2::element_blank(),axis.line.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                            text = ggplot2::element_text(size=30), legend.key = ggplot2::element_rect(colour = "white"), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, size = 25),
                            plot.subtitle = ggplot2::element_text(size = 15),axis.text.x = ggplot2::element_text(angle = 60,vjust = 1, hjust=1, size = 12), axis.text.y = ggplot2::element_text(size = 12)) +
                ggplot2::labs(title = title_out, x = "", y = "", fill = "") +
                ggplot2::scale_fill_viridis_c() +
                ggplot2::facet_wrap(~feature_name, scales = 'free')

    return(plot_out)
  }

  incidence_and_plot<- function(tree){
    g <- tree@tree
    incidence_df <- g %>% get_edge_incidence()
    if(method == 'overlap'){
      out_plot <- incidence_df %>% edge_overlap() %>% plot_overlap()

    }else if(method == 'jaccard'){
      out_plot <- incidence_df %>% edge_jaccard_similarity() %>% plot_overlap()

    }else{
      stop('Method not implemented yet!')
    }

    return(out_plot)
  }

  group.by <- get_feature_names(trees, group.by)

  if(inherits(trees, 'list')){
    out_plots <- list()
    for(i in 1:length(trees)){
      out_plots[[i]] <- lapply(trees[[i]], function(x)  x %>% incidence_and_plot())
    }

  }else if(inherits(trees, 'AntibodyForests')){
    out_plots <- trees %>% incidence_and_plot()

  }else{
    stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
  }

  return(out_plots)
}
