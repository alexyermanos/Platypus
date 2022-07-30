#' Plots the resulting node metrics from the AntibodyForests_metrics function

#'@description Will plot the resulting node metrics from the AntibodyForests_metrics function either as violin plots (plot.format = 'violin' or as a scatter plot of the principal components of the metrics dataframe - as long as multiple node metrics are calculated per node).
#' Requires the AntibodyForests_plot_metrics to be called before (as this function uses the resulting node_metrics dataframes).

#' @param trees nested list of AntibodyForests objects or single object, as obtained from the AntibodyForests function.
#' @param plot.format string - 'violin' for violin plots of node metrics per feature (as determined by the feature parameter), or 'pca' for performing PCA on the node metrics dataframe.
#' @param metrics.to.plot vector of strings - the metrics to be plotted from the metrics dataframe (must be already calculated using the AntibodyForests_metrics function)
#' @param group.by string - whether to group the violin/scatter plots by additional features (e.g., group.by = 'sample_id' to get violin plots per feature, per each unique sample).
#' @param max.groups integer - maximum number of groups to be considered in the resulting plots if group.by is not NULL.
#' @param specific.groups vector of strings - specific groups to plot if group.by is not NULL.
#' @param sample.by string - additional grouping factor (e.g., group.by can be set to 'clonotype_id' and sample.by to 'sample_id' for plots grouped by both clonotypes and samples)
#' @param features string - will determine the point colors in the PCA scatterplot.

#' @return either a violin plot or a scatter plot of the node metrics, as specified in the plot.format parameter
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_plot_metrics(trees,
#' plot.format = 'violin', metrics.to.plot = 'degree',
#' group.by = 'sample_id', sample.by = NULL)
#'}



AntibodyForests_plot_metrics <- function(trees,
                                         plot.format,
                                         metrics.to.plot,
                                         group.by,
                                         max.groups,
                                         specific.groups,
                                         sample.by,
                                         features
                                         ){

  if(missing(trees)) stop('Please input a nested list of AntibodyForests objects after having called ')
  if(missing(plot.format)) plot.format <- 'violin'
  if(missing(metrics.to.plot)) metrics.to.plot <- NULL
  if(missing(group.by)) group.by <- NULL
  if(missing(max.groups)) max.groups <- NULL
  if(missing(specific.groups)) specific.groups <- NULL
  if(missing(sample.by)) sample.by <- NULL
  if(missing(features)) features <- NULL


  node_metrics <- c('betweenness', 'closeness', 'eigenvector', 'authority_score',
                    'local_cluster_coefficient', 'average_cluster_coefficient',
                    'strength', 'degree',
                    'eccentricity', 'pagerank', 'daughters',
                    'path_from_germline', 'weighted_path_from_germline',
                    'path_from_most_expanded',
                    'weighted_path_from_most_expanded',
                    'path_from_hub', 'weighted_path_from_hub')


  if(is.null(metrics.to.plot)){
    pca_metrics <- c(node_metrics, 'distance_from_germline', 'cell_number')
  }else{
    pca_metrics <- metrics.to.plot
  }

  if(is.null(metrics.to.plot)){
    metrics.to.plot <- c('degree')
  }

  if(plot.format == 'pca'){
    group.by <- NULL
  }


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

  get_metrics_mean <- function(metrics_df){

    use.mean <- NULL
    plot.level = "sample" #setting defaults as now parameter is listed in funct def
    sample_id <- NULL
    clonotype_id <- NULL
    metric_q <- NULL
    median <- NULL

    if(!is.null(use.mean) & plot.level != 'clonotype'){

      for(metric in metrics.to.plot){
        if(plot.level == 'sample'){
          metrics_df <- metrics_df %>%
                        dplyr::group_by(sample_id, clonotype_id)

        }else if(plot.level == 'global'){
          metrics_df <- metrics_df %>%
                        dplyr::group_by(sample_id)
        }

        if(use.mean == 'avg'){
          metrics_df <- metrics_df %>% dplyr::mutate(metric_name = mean(!!metric_q))
          metrics_df[[metric]] <- NULL
        }else if(use.mean == 'median'){
          metrics_df <- metrics_df %>% dplyr::mutate(metric_name = median(!!metric_q))
          metrics_df[[metric]] <- NULL

        }
        names(metrics_df)[names(metrics_df) == 'metric_name'] <- metric
      }

      metrics_df <- metrics_df %>% dplyr::ungroup()
    }

    return(metrics_df)
  }

  get_metrics_df <- function(trees){

    cell_number <- NULL

    if(inherits(trees, 'list')){
      metrics_dfs <- vector(mode = 'list', length = length(trees))

      for(i in 1:length(trees)){
        if(inherits(trees[[i]], 'list')){
          dfs <- lapply(trees[[i]], function(x) x@metrics)
          metrics_dfs[[i]] <- do.call('rbind', dfs)
        }else{
          metrics_dfs[[i]] <- trees[[i]]@metrics
        }
      }

      metrics_df <- do.call('rbind', metrics_dfs)

    }else{
      metrics_df <- trees@metrics
    }


    metrics_df$name <- 1:nrow(metrics_df)
    metrics_df$label <- 1:nrow(metrics_df)

    if(is.null(group.by)){
      metrics_df$global <- 'global'
      group.by <- 'global'
    }
    if(!is.null(sample.by)){
      metrics_df <- split(metrics_df, metrics_df[sample.by])
    }else{
      metrics_df <- list(metrics_df)
      names(metrics_df) <- 'global'
    }

    if(!is.null(group.by)){
      if(!is.null(specific.groups)){
        for(i in 1:length(metrics_df)){
          metrics_df[[i]] <- metrics_df[[i]][metrics_df[[i]][[group.by]] %in% specific.groups,]

        }
      }

      if(!is.null(max.groups)){
        for(i in 1:length(metrics_df)){
          metrics_df[[i]] <- metrics_df[[i]] %>%
                             dplyr::group_by_at(group.by) %>%
                             dplyr::mutate(cells_per_group = sum(cell_number)) %>%
                             dplyr::ungroup()

          temp <- metrics_df[[i]] %>% dplyr::distinct_at(group.by, .keep_all = T)
          temp <- cbind(temp[group.by], temp$cells_per_group)
          temp <- temp[order(temp[,2], decreasing = T),]
          chosen_groups <- temp[1:max.groups,1]

          metrics_df[[i]] <- metrics_df[[i]][metrics_df[[i]][[group.by]] %in% chosen_groups,]
        }
      }
    }

    return(metrics_df)
  }


  bar_plot <- function(tree_subset, sample_id, clonotype_id, color.by){

      metrics_df <- tree_subset %>% get_metrics_df() %>% get_metrics_mean()

      final_plots <- vector(mode = 'list', length = length(metrics.to.plot))

      for(i in 1:length(metrics.to.plot)){
        final_plots[[i]] <- vector(mode = 'list', length = length(color.by))
        for(j in 1:length(color.by)){

          final_plots[[i]][[j]] <- ggplot2::ggplot(data = metrics_df, ggplot2::aes_string(x = 'clonotype_id', y = metrics.to.plot[i]), color = 'VDJ_cgene') +
              ggplot2::geom_bar(stat="identity", width=0.6, color="black") +
              #ggplot2::scale_y_continuous(expand = c(0,0)) +
              ggplot2::theme_bw() +
              ggplot2::theme_classic() +
              ggplot2::labs(col = color.by[j]) +
              ggplot2::labs(title = paste0(sample_id, ' - ', clonotype_id)) +
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x =  ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))

            }
        }

      #final_plots <- cowplot::plot_grid(plotlist = final_plots, align = 'h')

      return(final_plots)
  }

  pca_plot <- function(metrics_df, features, pca_metrics){

    PC1 <- NULL
    PC2 <- NULL
    feature <- NULL

    pca_matrix <- as.matrix(metrics_df[,which(colnames(metrics_df) %in% pca_metrics)])
    pca_matrix[is.na(pca_matrix)] <- 0
    pca_matrix[is.nan(pca_matrix)] <- 0

    pca <- stats::prcomp(t(pca_matrix), scale = F)
    df <- data.frame(PC1 = unname(pca$rotation[,1]), PC2 = unname(pca$rotation[,2]))

    final_plots <- vector(mode = 'list', length = length(features))

    for(i in 1:length(features)){
      plot_df <- df
      plot_df$feature <- metrics_df[,features[i]]
      plot_df$feature[is.na(plot_df$feature) | plot_df$feature==''] <- 'unknown'
      plot_df$feature_name <- features[i]

      final_plots[[i]] <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = PC1, y = PC2, color = stats::reorder(feature, nchar(feature)))) + #optional: size = cell_number
              ggplot2::geom_point(size = 0.75) +
              ggplot2::theme_bw() +
              ggplot2::theme_classic() +
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::labs(col = features[i], title = features[i])

    }
    #final_plots <- cowplot::plot_grid(plotlist = final_plots, align = 'h')
    names(final_plots) <- features
    return(final_plots)
  }


  violin_plot <- function(metrics_df, metrics.to.plot){

      final_plots <- vector(mode = 'list', length = length(metrics.to.plot))
      metrics_df[,group.by][metrics_df[,group.by] == '' | is.na(metrics_df[,group.by])] <- 'unknown'

      for(i in 1:length(metrics.to.plot)){

        final_plots[[i]] <- ggplot2::ggplot(data = metrics_df, ggplot2::aes_string(x = group.by, y = metrics.to.plot[i])) +
                            ggplot2::geom_violin(ggplot2::aes_string(fill = group.by), trim = F) +
                            ggplot2::geom_boxplot(width=0.1, fill = 'white') +
                            ggplot2::theme_bw() +
                            ggplot2::theme_classic() +
                            ggplot2::labs(title = metrics.to.plot[i]) +
                            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x =  ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
      }

      #final_plots <- cowplot::plot_grid(plotlist = final_plots, align = 'h')
      names(final_plots) <- metrics.to.plot
      return(final_plots)
  }

  metrics_dfs <- get_metrics_df(trees)
  output_plots <- vector(mode = 'list', length = length(metrics_dfs))

  features <- get_feature_names(trees, features)
  for(j in 1:length(metrics_dfs)){
    if(plot.format == 'violin'){
      output_plots[[j]] <- violin_plot(metrics_dfs[[j]], metrics.to.plot = metrics.to.plot)
    }else if(plot.format == 'pca'){
      in_df <- metrics_dfs[[j]]
      output_plots[[j]] <- pca_plot(metrics_df = in_df, features = features, pca_metrics = pca_metrics)
    }else{
      stop('Plot format is not available/recognized')
    }
  }

  #names(output_plots) <- names(metrics_dfs)
  return(output_plots)
}
