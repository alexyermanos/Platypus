#'This is a function to infer single cell trajectories and identifying lineage structures on clustered cells. Using the slingshot library
#'@param GEX GEX output of the VDJ_GEX_matrix function (VDJ_GEX_matrix[[2]]))
#'@param cluster.num A seurat cluster number for starting point of the lineage. Can be identified by using Seurat::DimPlot(VGM[[2]],group.by = "seurat_clusters"). Default is "0".
#'@param grouping Determine by which identifier to group by. E.g. 'group_id' or default 'seurat_clusters' which are automatically generated in the clustering process.
#'@return Returns a list. Element [[1]] returns updated GEX object with the inferred pseudotime trajectories per lineage. [[2]] returns the UMAP with the grouped cells. [[3]] and [[4]] show the slingshot inferred trajectories in two different styles.
#'@export
#'@examples
#' \donttest{
#' try({
#' lineage_trajectories <- GEX_lineage_trajectories(Platypus::small_vgm[[2]],
#'  grouping = 'group_id',
#'  cluster.num = "3")
#'})
#'}

GEX_lineage_trajectories <- function( GEX,
                                      grouping,
                                      cluster.num){

  if(missing(GEX)) stop("Please provide GEX input for this function")
  if (missing(grouping)){
    grouping <- 'seurat_clusters' # as default
  }
  if(missing(cluster.num)){ cluster.num <- 0}# as default

  platypus.version <- "v3"

  ########
  seu <- NULL
  sce <- NULL
  sds <- NULL
  cell_colors_clust <- NULL
  lineages <- NULL

  ########
  message("Generate UMAP plot...")
  umap_plot <- Seurat::DimPlot(GEX, reduction = "umap",
                               group.by = grouping, pt.size = 0.5, label = TRUE)

  #clustering for slinghsot
  seu <- Seurat::RunUMAP(GEX, dims = 1:50)


  message("Convert data into a SingleCellExperiment object for the next steps...")
  sce<- Seurat::as.SingleCellExperiment(GEX)


  message("Infer slingshot lineage trajectories...")
  sds <- slingshot::slingshot(sce,
                              clusterLabels = seu$seurat_clusters,
                              reducedDim = 'UMAP',
                              start.clus = cluster.num,
                              stretch = 0)

  #slingshot doesn't support ggplot. Create own colour palette
  cell_pal <- function(cell_vars, pal_fun,...) {
    if (is.numeric(cell_vars)) {
      pal <- pal_fun(100, ...)
      return(pal[cut(cell_vars, breaks = 100)])
    } else {
      categories <- sort(unique(cell_vars))
      pal <- stats::setNames(pal_fun(length(categories), ...), categories)
      return(pal[cell_vars])
    }
  }

  #color clusters
  cell_colors_clust <- cell_pal(seu$seurat_clusters, scales::hue_pal())


  message("Generate lineage trajectories plot...")

  invisible(plot(SingleCellExperiment::reducedDims(sds)$UMAP, col = cell_colors_clust, pch = 16, cex = 0.5))
  invisible(graphics::lines(slingshot::SlingshotDataSet(sds), lwd = 1, type = 'lineages', col = 'black'))
  lineage_plot= invisible(grDevices::recordPlot())

  invisible(plot(SingleCellExperiment::reducedDims(sds)$UMAP, col = cell_colors_clust, pch = 16, cex = 0.5))
  invisible(graphics::lines(slingshot::SlingshotDataSet(sds), lwd = 1,  col = 'black'))
  lineage_line_plot = invisible(grDevices::recordPlot())


  #update GEX with the inferred slingshot pseudotimes per lineages
  lineages <- slingshot::slingPseudotime(sds)
  for (i in 1:length(colnames(lineages))){
    GEX <- SeuratObject::AddMetaData(GEX,lineages[, i], col.name = paste0('pseudotime_lineage',i))
  }


  message("DONE")
  return(list(GEX, umap_plot, lineage_plot, lineage_line_plot))

}
