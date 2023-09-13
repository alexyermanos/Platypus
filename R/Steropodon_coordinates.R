#' Obtain a feature matrix of aligned coordinates across all Steropodon objects


#'@description Function used to obtain a feature matrix (n_structures x aligned_coordinates) for superposed Steropodon structures.
#' The feature matrix can be visualized as a 2D scatter plot, using the dimensionality reduction algorithm specified by dim.reduction, and colored by the VGM features mentioned in the color.by parameter.
#' Otherwise, it will output a 2D Matrix object for downstream analysis (e.g., clustering using Steropodon_cluster).

#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param use.pdbs bool - if TRUE, will use coordinates of the bio3d pdbs object obtained from Steropodon_sequence_alignment function (aligned coordinate vectors). If use.pdbs is set to FALSE, the resulting feature matrix will include all unaligned coordinate vectors, matching the length of the largest coordinate vector from the Steropodon nested list (with NA values added to match lengths).
#' @param feature.imputation string - method to impute NA values in the aligned coordinate matrix: 'zeros' will add zeros, 'mean' will get the mean across each sample, 'median' for the median value across each sample row, 'iterative_pca' will use the missMDA::imputePCA() function to infer the missing values.
#' @param remove.na bool - if TRUE, will remove NA values from the aligned coordinate matrix.
#' @param scale bool - if TRUE, will scale the output coordinate/feature matrix using Seurat::ScaleData().
#' @param variable.features integer or NULL - will select N most variable features in the coordinate matrix using Seurat::FindVariableFeatures(matrix, nfeatures = variable.features).
#' @param plot.results bool - if TRUE, will output a 2D scatter plot of the coordinate matrix after dimensionality reduction (defined using the dim.reduction parameter). Points can be colored by a custom VGM label specified in the color.by parameter. The VGM also needs to be included in the VGM parameter.
#' @param dim.reduction string - 'pca', 'tsne', or 'umap' for the dimensionality reduction algorithm performed before plotting.
#' @param color.by string or vector of strings - the column name in the VDJ/GEX parts of the VGM object which will be used to color each Steropodon structure/ scatter plot point.
#' @param VGM the VGM object to select the custom labels from if color.by is not NULL.
#' @param cluster.method string or NULL - the bluster clustering method for the 2D scatter plots. For more information, see the bluster::clusterRows() function.
#' @param additional.dim.reduction.parameters named list - additional parameters for the dimensionality reduction algorithms.
#' @param per.cell bool - if TRUE, will create the output scatter plot/ feature matrix per cell (barcode x coordinates) by pooling all cells assigned to a Steropodon-modelled receptor sequence.


#' @return a feature/coordinate matrix for downstream analyses or a 2D scatter plot.
#' @export
#' @examples
#' \dontrun{
#' superposed_w_pdbs <- superposed_seq_struct %>%
#' Steropodon_sequence_alignment(structure = 'structure')
#' superposed_w_pdbs %>% Steropodon_coordinates(use.pdbs = T,
#' color.by = c('site', 'age', 'VDJ_vgene'),
#' dim.reduction = 'pca', VGM = VGM)
#'}

Steropodon_coordinates <- function(steropodon.object,
                                   structure,
                                   use.pdbs,
                                   feature.imputation,
                                   remove.na,
                                   scale,
                                   variable.features,
                                   plot.results,
                                   dim.reduction,
                                   color.by,
                                   VGM,
                                   cluster.method,
                                   additional.dim.reduction.parameters,
                                   per.cell
                                  ){

  if(missing(steropodon.object)) stop('Please input your Steropodon object!')
  if(missing(structure)) structure <- 'structure'
  if(missing(use.pdbs)) use.pdbs <- T #recommended to get aligned coordinates as features (for dim reduction)
  if(missing(feature.imputation)) feature.imputation <- NULL
  if(missing(remove.na)) remove.na <- T
  if(missing(scale)) scale <- T
  if(missing(variable.features)) variable.features <- 100
  if(missing(plot.results)) plot.results <- T
  if(missing(dim.reduction)) dim.reduction <- 'tsne'
  if(missing(color.by)) color.by <- 'sample_id'
  if(missing(VGM)) VGM <- NULL
  if(missing(cluster.method)) cluster.method <- NULL #will cluster on dim reduced, not on whole features (as opposed to Steropodon_cluster)
  if(missing(additional.dim.reduction.parameters)) additional.dim.reduction.parameters <- list()
  if(missing(per.cell)) per.cell <- F

  get_feature_matrix <- function(steropodon.list,
                                 structure,
                                 use.pdbs,
                                 remove.na
                                 ){

    if(!use.pdbs){
      coords_vector <- lapply(steropodon.list, function(x) {struct <- select_structure(x, structure = structure)
                                                            coord <- c(struct$xyz)
                                                            return(t(as.matrix(coord)))
                                                           })

      coords <- plyr::rbind.fill.matrix(coords_vector)
      rownames(coords) <- names(steropodon.list)

    }else{
      pdbs <- steropodon.list[[1]]@pdbs
      if(is.null(pdbs)){
        stop('Could not find the aligned PDBs object. Please use Steropodon_sequence_alignment before!')
      }
      coords <- pdbs$xyz[]
    }

    if(remove.na){
      coords <- t(na.omit(t(coords)))
    }

    return(coords)
  }

  impute_feature_matrix <- function(feature.matrix,
                                    feature.imputation){

    if(feature.imputation == 'zeros'){
      feature.matrix[is.na(feature.matrix)] <- 0
    }else if(feature.imputation == 'mean'){
      feature.matrix <- missMethods::impute_mean(feature.matrix)

    }else if(feature.imputation == 'median'){
      feature.matrix <- missMethods::impute_median(feature.matrix)

    }else if(feature.imputation == 'iterative_pca'){
      nPCs <- missMDA::estim_ncpPCA(feature.matrix)
      feature.matrix <- missMDA::imputePCA(feature.matrix, method = 'EM', ncp = nPCs$ncp, scale = FALSE)

    }else{
      stop('More feature imputation method will be added in the next updates!')
    }

    return(feature.matrix)

  }


  features_dim_reduction <- function(feature.matrix,
                                     steropodon.list,
                                     dim.reduction,
                                     color.by,
                                     VGM,
                                     additional.dim.reduction.parameters){

    if(dim.reduction == 'pca'){
     out <- stats::prcomp(t(feature.matrix))$rotation[,1:2]

    }else if(dim.reduction == 'tsne'){
     params <- list(X = feature.matrix, check_duplicates = F)
     params <- c(params, additional.dim.reduction.parameters)
     out <- do.call(Rtsne::Rtsne, params)$Y
     rownames(out) <- rownames(feature.matrix)

    }else if(dim.reduction == 'umap'){
     params <- list(d = feature.matrix)
     params <- c(params, additional.dim.reduction.parameters)
     out <- do.call(umap::umap, params)$layout

    }else{
     stop('Dim reduction algorithm is not available!')
    }

    out <- as.data.frame(out)
    colnames(out) <- c('X', 'Y')

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

  out <- get_feature_matrix(steropodon.list = steropodon_list,
                            structure = structure,
                            use.pdbs = use.pdbs,
                            remove.na = remove.na)

  if(!remove.na & !is.null(feature.imputation)){
    out <- impute_feature_matrix(out,
                                 feature.imputation = feature.imputation)
  }

  if(per.cell){
    barcode_list <- lapply(steropodon_list, function(x) x@barcodes)
    barcode_counts <- lapply(barcode_list, function(x) length(x))
    out <- as.data.frame(out)
    out <- out[rep(1:nrow(out), times = unlist(barcode_counts)), ]
    out <- as.matrix(out)
    rownames(out) <- unname(unlist(barcode_list))
  }

  #Use the Seurat utilities for scaling/ finding variable features
  if(!is.null(variable.features)){
    out <- Seurat::CreateAssayObject(counts = t(out))
    out <- Seurat::FindVariableFeatures(out, nfeatures = variable.features)
    out <- t(as.matrix(Seurat::GetAssayData(out)))
  }

  if(scale){
    out <- Seurat::CreateAssayObject(counts = t(out))
    out <- Seurat::ScaleData(out)
    out <- t(as.matrix(Seurat::GetAssayData(out)))
  }

  if(plot.results & !per.cell){
    out <- features_dim_reduction(feature.matrix = out,
                                  steropodon.list = steropodon_list,
                                  dim.reduction = dim.reduction,
                                  color.by = color.by,
                                  VGM = VGM,
                                  additional.dim.reduction.parameters = additional.dim.reduction.parameters)

    out <- plot_embeddings(out)
  }


  return(out)
}
