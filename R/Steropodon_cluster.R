#' Structural clustering from Steropodon inferred coordinates/feature matrix or distance matrices


#'@description Function for clustering superposed coordinates from Steropodon_coordinates (feature matrices) or distance matrices from Steropodon_distance.
#' Different clustering algorithms can be used from the bluster package (https://bioconductor.org/packages/release/bioc/vignettes/bluster/inst/doc/clusterRows.html).
#' Specific clustering parameters (to finetune bluster's functions) can be specified as a named list in additional.cluster.parameters.
#' Dimensionality reduction (i.e., PCA, T-SNE, UMAP) can be performed before clustering, as defined in the dim.reduction parameter.
#' Additionally, plot.dim.reduction defines the dimensionality reduction algorithm for the 2D scatterplots, if plot.results is set to TRUE.

#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param distance.matrix a distance matrix obtained from Steropodon_distances for clustering.
#' @param feature.matrix a feature/superposed coordinate matrix obtained from Steropodon_coordinates for clustering. Ensure structures are initially aligned using Steropodon_superpose.
#' @param multimodal.object a multimodal Seurat object as obtained from Steropodon_multimodal (integrated GEX and structural features).
#' @param dim.reduction string - 'pca', 'tsne', or 'umap' for the dimensionality reduction algorithm performed before clustering.
#' @param plot.dim.reduction string - 'pca', 'tsne', or 'umap' for the dimensionality reduction algorithm performed after clustering/for plotting if plot.results is set to TRUE.
#' @param cluster.method bluster function - see the bluster package documentation (https://bioconductor.org/packages/release/bioc/vignettes/bluster/inst/doc/clusterRows.html) for all clustering algorithms available (for clustering either feature or distance matrices).
#' @param additional.cluster.parameters named list - additional parameters for the clustering algorithms.
#' @param additional.dim.reduction.parameters named list - additional parameters for the dimensionality reduction algorithms.
#' @param plot.results boolean - if TRUE, will output a 2D scatter plot for all structures in the Steropodon nested list, colored by cluster ID. Else, will output the Steropodon objects with a new 'cluster' slot/attribute.


#' @return a scatter plot colored by cluster ID or a nested list of Steropodon objects (input in the steropodon.object parameter) with a new 'cluster' slot.
#' @export
#' @examples
#' \dontrun{
#' feature_clusters <- Steropodon_cluster(steropodon.object = steropodon_trimmed,
#' feature.matrix = feature_matrix,
#' dim.reduction = 'tsne',
#' plot.dim.reduction = 'tsne',
#' cluster.method = bluster::HclustParam,
#' additional.dim.reduction.parameters = list(perplexity = 10),
#' plot.results = T)
#'}


Steropodon_cluster <- function(steropodon.object,
                               distance.matrix,
                               feature.matrix,
                               multimodal.object,
                               dim.reduction,
                               plot.dim.reduction,
                               cluster.method,
                               additional.cluster.parameters,
                               additional.dim.reduction.parameters,
                               plot.results
                              ){

   if(missing(steropodon.object)) stop('Please input your Steropodon object')
   if(missing(distance.matrix)) distance.matrix <- NULL
   if(missing(feature.matrix)) feature.matrix <- NULL
   if(missing(multimodal.object)) multimodal.object <- NULL
   if(missing(dim.reduction)) dim.reduction <- 'pca'
   if(missing(plot.dim.reduction)) plot.dim.reduction <- 'pca'
   if(missing(cluster.method)) cluster.method <- NULL
   if(missing(additional.cluster.parameters)) additional.cluster.parameters <- list()
   if(missing(additional.dim.reduction.parameters)) additional.dim.reduction.parameters <- list()
   if(missing(plot.results)) plot.results <- T

   X <- NULL
   Y <- NULL
   color.by <- NULL

   plot_embeddings <- function(embed.dfs){

     out_plots <- vector(mode = 'list', length = length(embed.dfs))

     for(i in 1:length(embed.dfs)){
       df <- embed.dfs[[i]]
       out_plots[[i]] <- ggplot2::ggplot(data = df, ggplot2::aes(x = X, y = Y, color = color.by)) +
                         ggplot2::geom_point(size = 1) +
                         ggplot2::theme_bw() +
                         ggplot2::theme_classic() +
                         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
                         ggplot2::labs(title = unique(df$feature.name), colour = unique(df$feature.name))
     }

     return(out_plots)
   }

   cluster_feature_matrix <- function(feature.matrix,
                                      dim.reduction,
                                      cluster.method,
                                      additional.cluster.parameters,
                                      additional.dim.reduction.parameters,
                                      plot.dim.reduction,
                                      plot.results){

      if(missing(feature.matrix)) stop('Please input your feature/coordinate matrix, as obtained from the Steropodon_coordinates function!')
      if(missing(dim.reduction)) dim.reduction <- 'pca'
      if(missing(cluster.method)) cluster.method <- bluster::HclustParam
      if(missing(additional.cluster.parameters)) additional.cluster.parameters <- NULL
      if(missing(additional.dim.reduction.parameters)) additional.dim.reduction.parameters <- list()
      if(missing(plot.dim.reduction)) plot.dim.reduction <- 'pca'
      if(missing(plot.results)) plot.results <- T

      if(!is.null(dim.reduction)){
        if(dim.reduction == 'pca'){
          feature.matrix <- stats::prcomp(t(feature.matrix))$rotation[,1:2]

        }else if(dim.reduction == 'tsne'){
          params <- list(X = feature.matrix, check_duplicates = F)
          params <- c(params, additional.dim.reduction.parameters)
          out <- do.call(Rtsne::Rtsne, params)$Y
          rownames(out) <- rownames(feature.matrix)
          feature.matrix <- out

        }else if(dim.reduction == 'umap'){
          params <- list(d = feature.matrix)
          params <- c(params, additional.dim.reduction.parameters)
          feature.matrix <- do.call(umap::umap, params)$layout
        }
      }

      clusters <- bluster::clusterRows(feature.matrix, cluster.method(additional.cluster.parameters))
      clusters <- as.vector(clusters)
      names(clusters) <- rownames(feature.matrix)


      if(plot.results){
        if(is.null(dim.reduction)){
          if(plot.dim.reduction == 'pca'){
            feature.matrix <- stats::prcomp(t(feature.matrix))$rotation[,1:2]

          }else if(plot.dim.reduction == 'tsne'){
            params <- list(X = feature.matrix, check_duplicates = F)
            params <- c(params, additional.dim.reduction.parameters)
            out <- do.call(Rtsne::Rtsne, params)$Y
            rownames(out) <- rownames(feature.matrix)
            feature.matrix <- out

          }else if(plot.dim.reduction == 'umap'){
            params <- list(d = feature.matrix)
            params <- c(params, additional.dim.reduction.parameters)
            feature.matrix <- do.call(umap::umap, params)$layout
          }
        }

        out <- as.data.frame(feature.matrix)
        colnames(out) <- c('X', 'Y')

        out$color.by <- as.factor(clusters)
        out$feature.name <- 'clusters'

        temp <- list()
        temp[[1]] <- out
        out <- temp

        p <- plot_embeddings(out)
        return(list(clusters = clusters, plot = p))
      }

      return(list(clusters = clusters, plot = NULL))
   }

   cluster_distance_matrix <- function(distance.matrix,
                                       cluster.method,
                                       plot.results,
                                       additional.cluster.parameters
                                       ){

      if(cluster.method == 'hclust'){
        params <- list(d = stats::as.dist(distance.matrix))
        clust <- do.call(stats::hclust, c(params, additional.cluster.parameters))

        if(is.null(additional.cluster.parameters$k)){
          additional.cluster.parameters$k <- 3
        }

        clusters <- stats::cutree(clust, k = additional.cluster.parameters$k)


        if(plot.results){
          unique_clusters <- length(unique(clusters))
          p <- stats::as.dendrogram(clust)
          p <- p %>%
               dendextend::color_branches(k = unique_clusters) %>%
               dendextend::color_labels(k = unique_clusters) %>%
               dendextend::set("branches_lwd", rep(2, length(clusters))) %>%
               dendextend::set("labels_cex", rep(0.75, length(clusters)))
        }

      }else if(cluster.method == 'medoids'){
        params <- list(x = stats::as.dist(distance.matrix))

        if(is.null(additional.cluster.parameters$k)){
          additional.cluster.parameters$k <- 3
        }

        p <- do.call(cluster::pam, c(params, additional.cluster.parameters))
        clusters <- p$clustering

      }else{
        stop('Method not implemented yet!')
      }

      if(plot.results){
        return(list(clusters = clusters, plot = p))
      }

      return(list(clusters = clusters, plot = NULL))
   }

   #TO DO:
   cluster_multimodal_object <- function(multimodal.object,
                                         cluster.method,
                                         additional.cluster.parameters,
                                         additional.dim.reduction.parameters,
                                         plot.results
                                        ){
   }


   if(is.null(cluster.method) & !is.null(distance.matrix)){
     cluster.method <- bluster::HclustParam
   }

   if(is.null(cluster.method) & !is.null(feature.matrix)){
     cluster.method <- 'hclust'
   }


   if(!is.null(feature.matrix)){
     if(length(additional.cluster.parameters) == 0){
       additional.cluster.parameters <- NULL
     }
     out <- cluster_feature_matrix(feature.matrix = feature.matrix,
                                   dim.reduction = dim.reduction,
                                   cluster.method = cluster.method,
                                   additional.cluster.parameters = additional.cluster.parameters,
                                   additional.dim.reduction.parameters = additional.dim.reduction.parameters,
                                   plot.dim.reduction = plot.dim.reduction,
                                   plot.results = plot.results)

   }else if(!is.null(distance.matrix)){
     out <- cluster_distance_matrix(distance.matrix = distance.matrix,
                                    cluster.method = cluster.method,
                                    plot.results = plot.results,
                                    additional.cluster.parameters = additional.cluster.parameters)

   }else{
     stop('Please input either a feature/coordinates matrix from Steropodon_coordinates or a distance matrix from Steropodon_distances!')
   }

   clusters <- out$clusters
   p <- out$plot


   steropodon_list <- unnest_steropodon(steropodon.object)
   obj_list <- steropodon_list[names(clusters)]

   for(i in 1:length(obj_list)){
     obj_list[[i]]@cluster <- unname(clusters[i])
   }

   steropodon_list[names(clusters)] <- obj_list
   steropodon_object <- nest_steropodon(steropodon_list)

   if(!is.null(p)){
     if(!inherits(p, 'list')){
       plot(p)
     }else{
       plot(p[[1]])
     }
   }

   return(steropodon_object)
}
