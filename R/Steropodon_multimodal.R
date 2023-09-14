#' Integrate structural and transcriptomic features


#' @description Function to integrate gene expression data and receptor structure coordinates/structural features into a single embedding.
#' This embedding can be further projected and visualized if plot.results is set to TRUE. For dimensionality reduction, we have implemented UMAP, T-SNE, and PCA.
#'
#' @param VGM the Platypus VGM object, with both VDJ and GEX information. The GEX data (as a Seurat object) will be used for multimodal integration with the structural features.
#' @param feature.matrix the coordinate feature matrix obtained from Steropodon_coordinates.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param n.variable.features vector of integers - number of variable feature to be selected for the coordinate features and GEX matrix, respectively.
#' @param pca.dims vector of integers - the number of principal components to be used for UMAP projection for coordinate and GEX features, respectively.
#' @param multimodal.integration string - method for the multimodal integration of GEX and structure feature matrices. Currently, only the weighted nearest neighbour approach from Seurat is implemented (multimodal.integration = 'wnn').
#' @param additional.integration.parameters named list - additional parameters for the multimodal integration algorithm.
#' @param additional.clustering.parameters named list - additional parameters for the clustering algorithm (see Seurat::FindClusters for more information).
#' @param plot.results bool - if TRUE, will plot a 2D scatter plot of the projected multimodal/integrated features.
#' @param color.by string or vector of strings - the column name in the VDJ/GEX parts of the VGM object which will be used to color each scatter plot point.

#' @return bar plots if plot.format = 'bar' or line plots if plot.format = 'line' of the properties specified in feature.
#' @export
#' @examples
#' \dontrun{
#' per_cell_features <- superposed_w_pdbs %>% \
#' Steropodon_coordinates(use.pdbs = T, plot.results = F, per.cell = T)
#' Steropodon_multimodal(VGM = VGM[[2]],
#'                      feature.matrix = per_cell_features,
#'                      multimodal.integration = 'wnn',
#'                      pca.dims = c(10, 20),
#'                      n.variable.features = c(100, 2000),
#'                      color.by = c('site', 'age', 'sample_id', 'seurat_clusters', 'VDJ_vgene'))
#'}


Steropodon_multimodal <- function(
                                  VGM,
                                  feature.matrix,
                                  structure,
                                  n.variable.features,
                                  pca.dims,
                                  multimodal.integration,
                                  additional.integration.parameters,
                                  additional.clustering.parameters,
                                  plot.results,
                                  color.by
                                ){

  if(missing(VGM)) stop('Please input a Seurat object or a VGM object for the multimodal analysis!')
  if(missing(feature.matrix)) stop('Please input a per-cell feature matrix, obtained from Steropodon_coordinates with per.cell = T!')
  if(missing(structure)) structure <- 'structure'
  if(missing(n.variable.features)) n.variable.features <- c(50, 2000)
  if(missing(pca.dims)) pca.dims <- c(10, 20)
  if(missing(multimodal.integration)) multimodal.integration <- 'wnn'
  if(missing(additional.integration.parameters)) additional.integration.parameters <- list()
  if(missing(additional.clustering.parameters)) additional.clustering.parameters <- list()
  if(missing(plot.results)) plot.results <- T
  if(missing(color.by)) color.by <- 'sample_id'


  create_multimodal_object <- function(feature.matrix,
                                       seurat.obj,
                                       n.variable.features,
                                       multimodal.integration,
                                       pca.dims,
                                       additional.integration.parameters,
                                       additional.clustering.parameters,
                                       structure
                                      ){

    struct.obj <- SeuratObject::CreateAssayObject(t(feature.matrix))
    cells <- rownames(feature.matrix)
    #Subset obj to match
    seurat.obj <- subset(seurat.obj, cells = cells)
    struct.obj <- subset(struct, cells = colnames(seurat.obj))
    seurat.obj[[structure]] <- struct

    Seurat::DefaultAssay(seurat.obj) <- structure
    seurat.obj <- suppressWarnings(seurat.obj %>%
                  Seurat::FindVariableFeatures(nfeatures = n.variable.features[1]) %>%
                  Seurat::ScaleData() %>%
                  Seurat::RunPCA(reduction.name = 'struct.pca'))


    Seurat::DefaultAssay(seurat.obj) <- 'RNA'
    seurat.obj <- suppressWarnings(seurat.obj %>%
                  Seurat::FindVariableFeatures(nfeatures = n.variable.features[2]) %>%
                  Seurat::ScaleData() %>%
                  Seurat::RunPCA(reduction.name = 'seurat.pca'))

    seurat.obj <- Seurat::RunUMAP(seurat.obj , reduction = 'struct.pca', dims = 1:pca.dims[1], assay = structure,
                  reduction.name = 'struct.umap', reduction.key = 'structUMAP_')

    seurat.obj <- Seurat::RunUMAP(seurat.obj , reduction = 'seurat.pca', dims = 1:pca.dims[2], assay = 'RNA',
                  reduction.name = 'seurat.umap', reduction.key = 'rnaUMAP_')


    if(multimodal.integration == 'wnn'){
      params <- list(object = seurat.obj,
                     reduction.list = list("struct.pca", "seurat.pca"),
                     dims.list = list(1:pca.dims[1], 1:pca.dims[2]),
                     modality.weight.name = "struct.weight",
                     weighted.nn.name = 'struct.nn',
                     knn.graph.name = 'struct.wknn',
                     snn.graph.name = 'struct.wsnn'
                     )

      params <- c(params, additional.integration.params)
      seurat.obj <- do.call(Seurat::FindMultiModalNeighbors, params)

      seurat.obj <- Seurat::RunUMAP(seurat.obj,
                                   nn.name = "struct.nn",
                                   reduction.name = "wnn.umap",
                                   reduction.key = "wnnUMAP_")

      params <- list(object = seurat.obj,
                     graph.name = 'wsnn',
                     verbose = False)

      params <- c(params, additional.clustering.params)
      seurat.obj <- do.call(Seurat::FindClusters, params)

    }else{
      stop('Multimodal integration method not implemented yet!')
    }

    return(seurat.obj)
  }


  plot_multimodal_embeddings <- function(multimodal.obj, color.by, multimodal.integration){
    if(multimodal.integration == 'wnn'){
      p1 <- Seurat::DimPlot(multimodal.obj, reduction = 'rna.umap', group.by = color.by, label = TRUE,
                  repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
      p2 <- Seurat::DimPlot(multimodal.obj, reduction = 'struct.umap', group.by = color.by, label = TRUE,
                  repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
      p3 <- Seurat::DimPlot(multimodal.obj, reduction = 'wnn.umap', group.by = color.by, label = TRUE,
                  repel = TRUE, label.size = 2.5)

      p4 <- Seurat::DimPlot(multimodal.obj, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5)

      return(list(compared = p1+p2+p3, clusters = p4))
    }else{
      stop('Multimodal integration method not implemented yet!')
    }
  }


  if(inherits(VGM, 'list')){
    seurat_obj <- VGM[[2]]
  }else{
    seurat_obj <- VGM
  }

  if(!inherits(seurat_obj, 'Seurat')){
    stop('Please input a Seurat object or a VGM object for the multimodal analysis!')
  }


  multimodal_obj <- create_multimodal_object(feature.matrix = feature.matrix,
                                             seurat.obj = seurat_obj,
                                             n.variable.features = n.variable.features,
                                             multimodal.integration = multimodal.integration,
                                             pca.dims = pca.dims,
                                             additional.integration.parameters = additional.integration.parameters,
                                             additional.clustering.parameters = additional.clustering.parameters,
                                             structure = structure)

  if(plot.results){
    out <- plot_multimodal_embeddings(multimodal.obj = multimodal_obj,
                                      color.by = color.by,
                                      multimodal.integration = multimodal.integration
                                      )
    return(out)
  }

  return(multimodal_obj)
}


#source('Steropodon_coordinates_DONE.R')
#VGM <- readRDS('../data/tnfr2_s1.rds')
#steropodon.object <- readRDS('./steropodon_annotated/2022-08-18 19:06:04_Steropodon_Object.rds')
#feature.matrix <- Steropodon_coordinates(steropodon.object, use.pdbs = F, remove.na = T, plot.results = F, per.cell = T)
#Steropodon_multimodal(VGM, feature.matrix)
