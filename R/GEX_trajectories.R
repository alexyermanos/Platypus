#' This is a function which infers trajectories along ordered cells on dimensionality reduced data. It projects trajectrories on a dim. red. plot such as Umap. This uses Monocle3 or Monocle2.
#' @param GEX GEX output of the VDJ_GEX_matrix function (VDJ_GEX_matrix[[2]]))
#' @param color.cells.by Column name in SummarizedExperiment::colData(GEX). To decide how the cells are colored in the output plot. E.g. color.cells.by = 'group_id' the cells will be colored based on their group_id.
#' @param reduction.method Which method to use for dimensionality reduction for monocle3. Supports "UMAP", "tSNE", "PCA" or "LSI". Default value is "UMAP".
#' @param monocle.version Version of monocle. Either monocle2 or monocle3. Default is monocle3.
#' @param genes Takes a vector of genes (e.g. genes = c('CD3E', 'CD4', 'CD8A', 'CD44')) to highlight the expression of these genes in UMAP and in the trajectory plot in monocle3. Default is NULL.
#' @param cluster.method Monocle3 gives two clustering options: Using the Leiden or the Louvain algo. Default is louvain.
#' @param label.cell.groups Whether to label cells in each group according to the most frequently occurring label(s) (as specified by color_cells_by) in the group. If false, plot_cells() simply adds a traditional color legend. Default is TRUE
#' @param label.groups.by.cluster Instead of labeling each cluster of cells, place each label once, at the centroid of all cells carrying that label. Default is TRUE
#' @param labels.per.group  How many labels to plot for each group of cells. Default is 1
#' @param group.label.size Font size to be used for cell group labels. Default is 1
#' @param ordering.cells.method In monocle2 you can choose between selecting genes with high dispersion across cells for ordering cells along a trajectory (= 'high.dispersion'). Or order cells based on genes which differ between clusters, uses an unsupervised procedure called "dpFeature" (= 'differ.genes'). Defalut is "differ.genes"
#' @return Returns a list.For monocle3: Element [[1]] returns a cell data set object with a new column for the UMAP clustering. This will be used for the GEX_pseudotime_trajectory_plot() function. [[2]] contains a plot of the clusters. [[3]] contains also a cluster plot but with the inferred trajectories. For monocle2: [[1]] cell data set object. [[2]] Trajetory plot with cells coloured based on their states (important to choose root state for pseudotime plot). [[3]] Trajectory plot based on color.cells.by
#' @export
#' @examples
#' \dontrun{
#'
#' trajectory_output <- GEX_trajectories(GEX = vgm[[2]],
#'  reduction.method = "UMAP",
#'  color.cells.by = "group_id",
#'  labels_per_group = 2,
#'  group_label_size = 3)
#'
#'  #visualizing gene expressions
#'  interesting_genes = c("Cxcr6", "Il7r")
#'  genes_trajectories <- GEX_trajectories(GEX = VGM$GEX,
#'   color.cells.by = "group_id",
#'    genes = interesting_genes)
#'
#' ##monocle2 ! DEPRECATED !
#' #trajectory_output <- GEX_trajectories(GEX = vgm[[2]],
#' #  monocle.version = "monocle2",
#' #  ordering.cells.method = "high.dispersion")
#'  }

GEX_trajectories <- function(GEX,
                             color.cells.by,
                             reduction.method,
                             cluster.method,
                             genes,
                             label.cell.groups,
                             label.groups.by.cluster,
                             labels.per.group,
                             group.label.size,
                             monocle.version,
                             ordering.cells.method){

  #missing inputs and setting default values
  if(missing(monocle.version)) { monocle.version <- "monocle3"}
  if(missing(GEX)) stop("Please provide GEX input for this function")
  if(missing(reduction.method)){
    reduction.method <- "UMAP" #Default
    }
  if(missing(genes)){
    genes <- NULL #Default
  }
  if(missing(cluster.method)){
    cluster.method <- "louvain"
  }

  if(missing(label.cell.groups)){
    label.cell.groups <- TRUE
  }
  if(missing(label.groups.by.cluster)){
    label.groups.by.cluster <- TRUE
  }
  if(missing(labels.per.group)){
    labels.per.group <- 1
  }
  if(missing(group.label.size)){
    group.label.size <- 1
  }

  platypus.version <- "v3"

  # Variable definitions
  cds <- NULL
  plot.clusters <- NULL
  learned.plot <- NULL
  trajectory.plot <- NULL
  mean_expression <- NULL
  dispersion_empirical <- NULL
  dispersion_fit <- NULL

  ##############

  if (monocle.version == "monocle3"){

    print("converting GEX to a cell data set object for monocle3...")
    cds <- SeuratWrappers::as.cell_data_set(GEX)

    #if genes are provided: check if the provided gene names actually exist in the cds object
    if (!(is.null(genes))){
      for (i in 1:length(genes)){
        if(genes[i] %in% rownames(SummarizedExperiment::rowData(cds)) == FALSE){
          stop(paste("The gene provided in the 'genes' param: ", genes[i], ",doesn't exist in the input object"))
        }
      }
      #make the row names in cds recognizable for Monocle3
      SummarizedExperiment::rowData(cds)$gene_name <- rownames(cds)
      SummarizedExperiment::rowData(cds)$gene_short_name <- SummarizedExperiment::rowData(cds)$gene_name
    }


    print("plot the clustered data...")
    cds <- monocle3::cluster_cells(cds, cluster_method = cluster.method, reduction_method = reduction.method) #grouping cells (using Leiden algo = community detection)

    if(missing(color.cells.by)){color.cells.by <- "seurat_clusters"}
    plot.clusters <- monocle3::plot_cells(cds,
                                          color_cells_by = color.cells.by,
                                          genes = genes,
                                          label_cell_groups= label.cell.groups ,
                                          label_groups_by_cluster = label.groups.by.cluster,
                                          labels_per_group = labels.per.group,
                                          group_label_size = group.label.size) # plotting the clusters



    #constructing single-cell trajectories
    print("infering trajectories by using monocle3. This step might take some while..")
    cds <- monocle3::learn_graph(cds)

    print("plot the trajectory plot...")
    learned.plot <- monocle3::plot_cells(cds,
                                         genes = genes,
                                         label_principal_points = TRUE,
                                         show_trajectory_graph=TRUE,
                                         label_cell_groups= TRUE ,
                                         label_groups_by_cluster = FALSE,)

    print("DONE")
    return(list(cds,plot.clusters, learned.plot))




  }else if (monocle.version == "monocle2"){

    stop("Monocle2 is deprecated. Please use monocle3")

    #if(missing(ordering.cells.method)){ordering.cells.method <- "differ.genes"}

    #print("converting GEX into a CellDataSet object")
    #cds <- Seurat::as.CellDataSet(GEX)
      #Size factor and dispersion will be used to normalize data and select genes for clustering
    #cds <- BiocGenerics::estimateSizeFactors(cds)
    #cds <- BiocGenerics::estimateDispersions(cds)

    #if(ordering.cells.method == "high.dispersion" ){
    #  disp_table <- monocle::dispersionTable(cds)
    #  clustering_genes <- subset(disp_table, mean_expression >= 0.5& dispersion_empirical >= 1 * dispersion_fit)$gene_id
    #  cds <- monocle::setOrderingFilter(cds, clustering_genes)
    #}
    #else if (ordering.cells.method == "differ.genes"){
    #  cds <- monocle::detectGenes(cds, min_expr = 0.1)
    #  Biobase::fData(cds)$use_for_ordering <-
    #    Biobase::fData(cds)$num_cells_expressed > 0.05 * ncol(cds)
    #}


    #print("reduce dimensions and cluster cells...")
    #cds <- monocle::reduceDimension(cds,
    #                       max_components = 2,
    #                       norm_method = 'log',
    #                       num_dim = 3,
    #                       reduction_method = 'tSNE',
    #                       verbose = T)

    #cds <- monocle::clusterCells(cds, verbose = F)

    #print("print clustered cells")
    #plot.clusters <- monocle::plot_cell_clusters(cds, color_by = color.cells.by)

    #print("Order cells and inferre trajectory. This step might take a while... ")

      #this function only worked for me when I loaded the DDRTree or monocle library before using the function
    #cds <- monocle::reduceDimension(cds, reduction_method = 'DDRTree')

    #cds <- monocle::orderCells(cds)

    #print("print trajectory plot")
    #if(missing(color.cells.by)){color.cells.by <- "State"}
    #trajectory.plot.state <- monocle::plot_cell_trajectory(cds, color_by = "State")
    #trajectory.plot <- monocle::plot_cell_trajectory(cds, color_by = color.cells.by)

    #print("DONE")
    #return(list(cds,plot.clusters, trajectory.plot.state, trajectory.plot))

  }


}










