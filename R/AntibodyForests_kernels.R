#' Graph kernel methods for graph structure/topology comparisons

#'@description Performs graph structural comparisons using graph kernel-based method. Currently available kernel methods include: the Weisfeiler-Lehman kernel, the graphlet kernel, and the random walk kernel.
#' @param trees (nested) list of AntibodyForests objects, as obtained from the AntibodyForests function, to be compared.
#' @param graph.type string - 'tree' will use the usual output of the AntibodyForests function (tree graphs), 'heterogeneous' will use the output of the AntibodyForests_heterogeneous function (bipartite networks for both cells and sequences).
#' @param kernel.method string - kernel method to be used, as implemented in the 'graphkernels' R package. 'weisfeiler_lehman' for the Weisfeiler-Lehman kernel, 'graphlet, and 'random_walk'.
#' @param additional.param integer - additional kernel options/parameters (e.g., kernel iterations for the Weisfeiler-Lehman kernel).
#' @param max.networks integer - maximum number of networks to be compared (will pick the networks with the most number of cells first).
#' @return Heatmap of the graph kernel values.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_kernels(trees, graph.type = 'tree', kernel.method = 'weisfeiler_lehman', additional.params = 10, max.networks = 50)
#'}


AntibodyForests_kernels <- function(trees,
                                    graph.type,
                                    kernel.method,
                                    additional.param,
                                    max.networks){

  if(missing(trees)) stop('Please input a list of AntibodyForests objects to be compared')
  if(missing(graph.type)) graph.type <- 'tree'
  if(missing(kernel.method)) kernel.method <- 'weisfeiler_lehman'
  if(missing(additional.param)) additional.param <- NULL
  if(missing(max.networks)) max.networks <- NULL

  graph_kernels <- function(tree.list){
    requireNamespace('graphkernels')
    requireNamespace('viridis')
    tree.list <- tree.list[order(nchar(names(tree.list)), names(tree.list))]

    if(!is.null(max.networks)){
      if(max.networks < length(tree.list)){
        tree.list <- tree.list[1:max.networks]
      }
    }

    for(i in 1:length(tree.list)){
      if(graph.type == 'tree'){
        g <- tree.list[[i]]@tree
      }else if(graph.type == 'heterogeneous'){
        g <- tree.list[[i]]@heterogeneous
      }
      igraph::V(g)$label <- 1:length(igraph::V(g))
      tree.list[[i]] <- g
    }


    if(kernel.method == 'weisfeiler_lehman'){
      if(is.null(additional.param)){
        additional.param <- 100
      }
      kernel_matrix <- suppressWarnings(graphkernels::CalculateWLKernel(tree.list, par=additional.param))
    }else if(kernel.method == 'graphlet'){
      if(is.null(additional.param)){
        additional.param <- 10
      }
      kernel_matrix <- suppressWarnings(graphkernels::CalculateGraphletKernel(tree.list, par=additional.param))
    }else if(kernel.method == 'random_walk'){
      if(is.null(additional.param)){
        additional.param <- 5
      }
      kernel_matrix <- suppressWarnings(graphkernels::CalculateKStepRandomWalkKernel(tree.list, par=additional.param))
    }else{
      stop('Kernel method not recognized!')
    }

    colnames(kernel_matrix) <- names(tree.list)
    rownames(kernel_matrix) <- names(tree.list)
    diag(kernel_matrix) <- 0
    #diag(kernel_matrix) <- max(kernel_matrix)
    kernel_heatmap <- pheatmap::pheatmap(kernel_matrix, color = viridis::viridis(10), cluster_rows=F, cluster_cols=F, main=paste0(kernel.method, ' kernel similarities'))


    return(kernel_heatmap)
  }

  if(inherits(trees[[1]], 'list')){
    plots <- vector(mode = 'list', length = length(trees))
    for(i in 1:length(trees)){
      plots[[i]] <- graph_kernels(trees[[i]])
    }

    names(plots) <- names(trees)
  }else{
    plots <- graph_kernels(trees)
  }

  return(plots)
}
