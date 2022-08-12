#' Bipartite sequence-cell networks in AntibodyForests

#'@description Creates a bipartite network from a Seurat object and an already inferred AntibodyForests sequence similarity/ minimum spanning tree network.
#' @param trees AntibodyForests object/list of AntibodyForests objects - the resulting sequence similarity or minimum spanning tree networks from the AntibodyForests function.
#' @param GEX.object Seurat object/ VGM[[2]] for the inferred AntibodyForests networks (must include all cells available in the AntibodyForests object).
#' @param node.features vector of strings - gene names in the Seurat object to be added as igraph vertex attributes in the resulting heterogeneous networks (will add gene expression per gene).
#' @param cell.graph.type string - graph algorithm for building the cell-cell/ nearest-neighbour graphs, as done by the Seurat::FindNeighbors() function. 'knn' for the K nearest-neighbour graphs, 'snn' for the shared nearest-neighbour graphs.
#' @param recluster.cells boolean - whether to recluster the resulting cell graphs or keep the already existing Seurat cluster definitions.
#' @param recluster.resolution numeric - recluster resolution for the Louvain algorithm if recluster.cells is TRUE.
#' @param snn.threshold numeric - SNN edge weight threshold to further prune edges in the cell graphs (increased value = sparser cell graphs). Defaults to 1/15 (as done in the Seurat::FindNeighbors function).
#' @param keep.largest.cc boolean - whether to keep only the largest connected component in the cell graphs or keep all singletons/doubletons/etc., as well.
#' @param parallel boolean - whether to execute the heterogeneous graph building algorithm in parallel or not. Requires the 'parallel' R package.

#' @return nested list of AntibodyForests objects for each clonotype and each sample/timepoint or a single AntibodyForests object, with a new added object slot for the heterogeneous graph.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_heterogeneous(trees, GEX.object = VGM[[2]], cell.graph.type = 'snn')
#'}

#Further analysis for bipartite - metacells on the cell graph, overlaps/all that in the resulting incidence matrix (now we ensure some overlap btw sequences and their corresponding cells)
#add incidence matrix analyses for metacell bipartite graphs

AntibodyForests_heterogeneous <- function(trees,
                                          GEX.object,
                                          node.features,
                                          cell.graph.type,
                                          recluster.cells,
                                          recluster.resolution,
                                          snn.threshold,
                                          keep.largest.cc,
                                          parallel){


  if(missing(trees)) stop('Please input your AntibodyForests object or a nested list of AntibodyForests objects!')
  if(missing(GEX.object)) stop('Please input your Seurat/GEX object')
  if(missing(node.features)) node.features  <- NULL
  if(missing(cell.graph.type)) cell.graph.type <- 'snn'
  if(missing(recluster.cells)) recluster.cells <- F
  if(missing(recluster.resolution)) recluster.resolution <- 0.5
  if(missing(snn.threshold)) snn.threshold <- 1/15
  if(missing(keep.largest.cc)) keep.largest.cc <- F
  if(missing(parallel)) parallel <- F



  create_cell_graph <- function(tree,
                                GEX_object = GEX.object,
                                node_features = node.features,
                                recluster_cells = recluster.cells,
                                recluster_resolution = recluster.resolution,
                                snn_threshold = snn.threshold,
                                cell_graph_type = cell.graph.type,
                                keep_largest_cc = keep.largest.cc){


    g <- tree@tree
    barcodes <- unlist(tree@barcodes)

    seurat_subset <- subset(GEX_object, cells = barcodes)
    seurat_subset <- Seurat::FindNeighbors(seurat_subset, prune.SNN = snn_threshold)

    if(recluster_cells){
      seurat_subset <- Seurat::FindClusters(seurat_subset, resolution = recluster_resolution)
    }

    if(cell_graph_type == 'knn'){
      adj_matrix <- seurat_subset@graphs$RNA_nn
    }else{
      adj_matrix <- seurat_subset@graphs$RNA_snn
    }

    adj_matrix <- as.matrix(adj_matrix)
    barcodes <- colnames(adj_matrix)
    colnames(adj_matrix) <- NULL
    rownames(adj_matrix) <- NULL
    diag(adj_matrix) <- 0

    clusters_per_cell <- unname(seurat_subset$seurat_clusters[barcodes])

    if(cell_graph_type == 'snn'){
      cell_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, weighted = T)
    }else{
      cell_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, weighted = NULL)
    }

    igraph::V(cell_graph)$cell_barcodes <- unlist(barcodes)
    igraph::V(cell_graph)$seurat_clusters <- unlist(clusters_per_cell)

    if(!is.null(node_features)){

      features <- Seurat::GetAssayData(seurat_subset)[node_features,]
      features <- as.data.frame(t(as.matrix(features)))

      for(feature in node_features){

        cell_graph <- igraph::set_vertex_attr(cell_graph, name = feature, value = unname(unlist(features[feature])))

      }
    }

    cell_graph <- igraph::delete_vertices(cell_graph, v = which(igraph::degree(cell_graph) == 0))

    if(!is.null(snn_threshold) & cell_graph_type == 'snn' & keep_largest_cc){
      ccs <- igraph::components(cell_graph)
      max_cc <- which.max(ccs$csize)
      cells_to_remove <- which(ccs$membership != max_cc)

      cell_graph <- igraph::delete_vertices(cell_graph, v = cells_to_remove)

    }

    return(list(sequence_graph = g, cell_graph = cell_graph))
  }
  create_heterogeneous_graph <- function(graph.list){

    cell_graph <- graph.list$cell_graph
    sequence_graph <- graph.list$sequence_graph

    cell_vertices <- igraph::as_data_frame(cell_graph, what = 'vertices')
    sequence_vertices <- igraph::as_data_frame(sequence_graph, what = 'vertices')

    cell_edges <- igraph::as_data_frame(cell_graph, what = 'edges')
    sequence_edges <- igraph::as_data_frame(sequence_graph, what = 'edges')

    cell_edges$from <- paste0('cell', cell_edges$from)
    cell_edges$to <- paste0('cell', cell_edges$to)
    cell_edges$type <- 'cell'

    sequence_edges$from <- paste0('sequence', sequence_edges$from)
    sequence_edges$to <- paste0('sequence', sequence_edges$to)
    sequence_edges$type <- 'sequence'

    cell_vertices$label <- paste0('cell', 1:nrow(cell_vertices))
    cell_vertices$type <- 'cell'
    cell_vertices$node_type <- 'cell'

    last_label <- max(sequence_vertices$label, na.rm = T)
    sequence_vertices$label[is.na(sequence_vertices$label)] <- (last_label+1):nrow(sequence_vertices)
    sequence_vertices$label <- paste0('sequence', sequence_vertices$label)
    sequence_vertices$type <- 'sequence'

    if(is.null(sequence_edges$weight)){
      sequence_edges$weight <- 1
    }

    heterogeneous_edges <- rbind(sequence_edges, cell_edges)

    sequence_barcodes <- sequence_vertices$cell_barcodes[!sequence_vertices$node_type %in% c('bulk', 'intermediate', 'germline', 'inferred')]
    names(sequence_barcodes) <- sequence_vertices$label[!sequence_vertices$node_type %in% c('bulk', 'intermediate', 'germline', 'inferred')]

    sequence_barcode_dfs <- list()
    for(i in 1:length(sequence_barcodes)){
      sequence_barcode_dfs[[i]] <- data.frame(label = rep(names(sequence_barcodes[i]), length(sequence_barcodes[[i]])),
                                              cell_barcodes = unlist(sequence_barcodes[[i]])
                                              )
    }


    sequence_barcode_df <- do.call('rbind', sequence_barcode_dfs)
    cell_barcode_df <- cell_vertices[names(cell_vertices) %in% c('cell_barcodes', 'label')]
    interlevel_edges <- merge(cell_barcode_df, sequence_barcode_df, by = 'cell_barcodes', all.x = T)
    interlevel_edges$cell_barcodes <- NULL
    interlevel_edges <- data.frame(from = c(interlevel_edges$label.x, interlevel_edges$label.y),
                                   to = c(interlevel_edges$label.y, interlevel_edges$label.x),
                                   weight = 1,
                                   type = 'interlevel')

    heterogeneous_edges <- rbind(heterogeneous_edges, interlevel_edges)


    sequence_cols <- unique(colnames(sequence_vertices))
    cells_cols <- unique(colnames(cell_vertices))

    missing_cols_sequences <- setdiff(cells_cols, sequence_cols)
    missing_cols_cells <- setdiff(sequence_cols, cells_cols)

    for(col in missing_cols_sequences){
      sequence_vertices[col] <- rep(NA, nrow(sequence_vertices))
    }

    for(col in missing_cols_cells){
      cell_vertices[col] <- rep(NA, nrow(cell_vertices))
    }

    heterogeneous_vertices <- rbind(sequence_vertices, cell_vertices)
    heterogeneous_vertices$name <- heterogeneous_vertices$label
    heterogeneous_graph <- igraph::graph_from_data_frame(d = heterogeneous_edges, directed = F, vertices = heterogeneous_vertices)
    #igraph::V(heterogeneous_graph)$seurat_clusters <- as.character(igraph::V(heterogeneous_graph)$seurat_clusters)
    igraph::V(heterogeneous_graph)$lvl <- ifelse(igraph::V(heterogeneous_graph)$type == 'sequence',1,2)
    return(heterogeneous_graph)
  }

  if(inherits(trees, 'list')){
    for(i in 1:length(trees)){
      heterogeneous_list <- vector(mode = "list", length = length(trees[[i]]))

      if(parallel){
        #requireNamespace('parallel')
        cores <- parallel::detectCores()
        heterogeneous_list <- parallel::mclapply(trees[[i]], mc.cores = cores,
                                                 FUN = function(x) {x %>% create_cell_graph() %>%  create_heterogeneous_graph()

                                                                    })

      }else{
        heterogeneous_list <- lapply(trees[[i]], function(x) {x %>% create_cell_graph() %>%  create_heterogeneous_graph()

                                                              })
      }
      for(j in 1:length(trees[[i]])){
        trees[[i]][[j]]@heterogeneous <- heterogeneous_list[[j]]
      }

    }

  }else if(inherits(trees, 'AntibodyForests')){
    trees@heterogeneous <- trees %>%
                           create_cell_graph() %>%
                           create_heterogeneous_graph()

  }else{
    stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
  }

  return(trees)
}
