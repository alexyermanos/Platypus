#DONE - some errors to fix/ find better WNN default params
Steropodon_multimodal <- function(
                                  VGM,
                                  feature.matrix,
                                  structural.name,
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
  if(missing(structural.name)) structural.name <- 'structure'
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
                                       structural.name
                                      ){

    struct.obj <- SeuratObject::CreateAssayObject(t(feature.matrix))
    cells <- rownames(feature.matrix)
    #Subset obj to match
    seurat.obj <- subset(seurat.obj, cells = cells)
    struct.obj <- subset(struct, cells = colnames(seurat.obj))
    seurat.obj[[structural.name]] <- struct

    Seurat::DefaultAssay(seurat.obj) <- structural.name
    seurat.obj <- suppressWarnings(seurat.obj %>%
                  Seurat::FindVariableFeatures(nfeatures = n.variable.features[1]) %>%
                  Seurat::ScaleData() %>%
                  Seurat::RunPCA(reduction.name = 'struct.pca'))


    Seurat::DefaultAssay(seurat.obj) <- 'RNA'
    seurat.obj <- suppressWarnings(seurat.obj %>%
                  Seurat::FindVariableFeatures(nfeatures = n.variable.features[2]) %>%
                  Seurat::ScaleData() %>%
                  Seurat::RunPCA(reduction.name = 'seurat.pca'))

    seurat.obj <- Seurat::RunUMAP(seurat.obj , reduction = 'struct.pca', dims = 1:pca.dims[1], assay = structural.name,
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
                                             n.variable.features = ,
                                             multimodal.integration = multimodal.integration,
                                             pca.dims = pca.dims,
                                             additional.integration.parameters = additional.integration.parameters,
                                             additional.clustering.parameters = additional.clustering.parameters,
                                             structural.name = structural.name)

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
