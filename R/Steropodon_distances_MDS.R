#' Performs multidimensional scaling on a distance matrix obtained from Steropodon_distances


#' @description Function for multidimensional scaling of the distance matrix obtained from Steropodon_distances.
#' Will also output a scatter plot colored by the VGM label column specified in the color.by parameter, if plot.results is set to TRUE.


#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param distance.matrix the distance matrix obtained from Steropodon_distances.
#' @param color.by string or vector of strings - the column name in the VDJ/GEX parts of the VGM object which will be used to color each Steropodon structure/ scatter plot point.
#' @param VGM  the VGM object to select the custom labels from if color.by is not NULL.
#' @param cluster.method string or NULL - the bluster clustering method for the 2D scatter plots. For more information, see the bluster::clusterRows() function.
#' @param plot.results bool - if TRUE, will output a 2D scatter plot of the distance matrix after MDS. Points can be colored by a custom VGM label specified in the color.by parameter. The VGM also needs to be included in the VGM parameter.


#' @return the MDS results or a 2D scatter plot if plot.results is set to TRUE.
#' @export
#' @examples
#' \dontrun{
#' distance_matrix <- superposed_seq_struct %>%
#' Steropodon_distances(distance.metric = 'rmsd', plot.results = F)

#' Steropodon_distances_MDS(steropodon.object = superposed_seq_struct,
#'                          distance.matrix = distance_matrix,
#'                          color.by = c('site', 'age', 'sample', 'VDJ_vgene'),
#'                          VGM = VGM)
#'}


Steropodon_distances_MDS <- function(steropodon.object,
                                     distance.matrix,
                                     color.by,
                                     VGM,
                                     cluster.method,
                                     plot.results
                                     ){

  if(missing(plot.results)) plot.results <- F
  if(missing(steropodon.object) & plot.results) stop('Please input your Steropodon object!')
  if(missing(steropodon.object) & !plot.results) steropodon.object <- NULL
  if(missing(distance.matrix)) stop('Please input the distance matrix obtained from Steropodon distances!')
  if(missing(color.by)) color.by <- 'sample'
  if(missing(VGM) & !plot.results) VGM <- NULL
  if(missing(VGM) & plot.results) stop('Please input your VGM object or VDJ matrix!')
  if(missing(cluster.method)) cluster.method <- NULL

  X <- NULL
  Y <- NULL

  distance_mds <- function(distance.matrix,
                            steropodon.list,
                            color.by,
                            VGM,
                            cluster.method,
                            plot.results
                            ){

   fit <- stats::cmdscale(distance.matrix, eig = TRUE, k = 2)

   if(!plot.results){
     return(fit$points)
   }

   out <- as.data.frame(fit$points)
   colnames(out) <- c('X', 'Y')

   steropodon.list <- steropodon.list[rownames(fit)]
   barcode_list <- lapply(steropodon.list, function(x) x@barcodes)

   if(!is.null(VGM)){
     if(inherits(VGM, 'list')){
       VGM <- VGM[[1]]
     }
   }

   if(is.null(color.by)) color.by <- 'sample_clonotype'

   out_dfs <- vector(mode = 'list', length = length(color.by))
   for(i in 1:length(color.by)){
     col <- color.by[i]
     if(col != 'sample' & col != 'sample_clonotype'){
       max_features <- lapply(barcode_list, function(x) {  features <- VGM[VGM$barcode %in% unlist(x),][[col]]
                                                           features[is.na(features) | features == ''] <- 'unknown'
                                                           feature_counts <- table(features)
                                                           max_feature <- names(feature_counts)[which.max(feature_counts)]
                                                           return(max_feature)
                                                         })
       out$color.by <- unlist(max_features)
       out$feature.name <- col

     }else{
       if(col == 'sample'){
         features <- lapply(rownames(out), function(x) stringr::str_split(x, '\\.')[[1]][1] )

       }

       if(col == 'sample_clonotype'){
         features <- lapply(rownames(out), function(x) paste0(stringr::str_split(x, '\\.')[[1]][1:2], collapse = '.'))

       }
       out$color.by <- unlist(features)
       out$feature.name <- col
     }

     out_dfs[[i]] <- out
   }

   if(!is.null(cluster.method)){
     out$color.by <- NULL
     out$feature.name <- NULL
     clusters <- bluster::clusterRows(out, cluster.method)
     out$color.by <- clusters
     out$feature.name <- 'clusters'
     out_dfs[[length(out_dfs) + 1]] <- out
   }

   return(out_dfs)
  }


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

  steropodon_list <- unnest_steropodon(steropodon.object)

  out <- distance_mds(distance.matrix = distance.matrix,
                      steropodon.list = steropodon_list,
                      color.by = color.by,
                      VGM = VGM,
                      cluster.method = cluster.method,
                      plot.results = plot.results
                     )

  if(plot.results){
    out <- plot_embeddings(out)
  }

  return(out)
}
