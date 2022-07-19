#'Converts the igraph networks of a given AntibodyForests object into a given (useful to convert the minimum spanning trees into a phylogenetic tree)

#'@description Will automatically convert the minimum spanning trees in a given AntibodyForests object into a phylogenetic tree as a phylo object. This new object will be added into the phylo slot of the AntibodyForests object.

#' @param trees nested list of AntibodyForests objects or single object, as obtained from the AntibodyForests function.
#' @param output.format string - 'treedata' will output the phylogenetic tree as a tidytree treedata object, 'phylo' as an ape::phylo object.
#' @param solve.multichotomies boolean - whether to remove multichotomies in the resulting phylogenetic tree using ape::multi2di
#' @param parallel boolean - whether to execute the main subroutine in parallel or not. Requires the 'parallel' R package to be installed.

#' @return nested list of AntibodyForests objects for each clonotype and each sample/timepoint or a single object, with a new phylo slot for the phylogenetic tree.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_phylo(trees, output.format = 'phylo')
#'}


AntibodyForests_phylo <- function(trees,
                                  output.format,
                                  solve.multichotomies,
                                  parallel){

  if(missing(trees)) stop('Please input your AntibodyForests object or a nested list of AntibodyForests objects!')
  if(missing(output.format)) output.format <- 'phylo'
  if(missing(solve.multichotomies)) solve.multichotomies <- F
  if(missing(parallel)) parallel <- T

  #requireNamespace('aphylo')
  #Function taken from the aphylo package as it could not be loaded properly
  phylo_object <- function(edge, tip.label, Nnode, edge.length = NULL, node.label = NULL, root.edge = NULL){
    structure(
      c(
        list(edge = edge),
        if (length(edge.length))
          list(edge.length = edge.length)
        else
          NULL,
        list(
          tip.label  = tip.label,
          Nnode      = Nnode,
          node.label = node.label
        ),
        if (length(root.edge))
          list(root.edge = root.edge)
        else
          NULL
      ),
      class = "phylo"
    )
  }

  to_phylo <- function(x, edge.length = NULL, root.edge = NULL){


    if (!inherits(x, "integer")) {
      label <- sort(unique(as.vector(x)))
      x <- matrix(match(as.vector(x), label), ncol=2L)

    }else{
      label <- range(as.vector(x))
      label <- label[1]:label[2]
    }

    nodes <- sort(unique(as.vector(x)))
    ideg <- tabulate(x[,2] - nodes[1L] + 1L, nbins = nodes[length(nodes)] - nodes[1L] + 1)
    odeg <- tabulate(x[,1] - nodes[1L] + 1L, nbins = nodes[length(nodes)] - nodes[1L] + 1)

    roots <- nodes[ideg == 0 & odeg > 0]
    leaves <- nodes[ideg == 1 & odeg == 0]
    inner <- nodes[ideg == 1 & odeg > 0]

    nodes <- c(leaves, roots, inner)

    iroots <- which(x[] == roots)
    lroots <- match(roots, nodes)

    ileaves <- which(x[] %in% leaves)
    lleaves <- match(x[ileaves], nodes)

    iinner <- which(x[] %in% inner)
    linner <- match(x[iinner], nodes)

    x[iroots] <- lroots
    x[ileaves] <- lleaves
    x[iinner] <- linner

    phylo_object(
      edge        = unname(x),
      edge.length = unname(edge.length),
      tip.label   = unname(label[leaves]),
      Nnode       = length(inner) + 1L,
      node.label  = unname(label[c(roots, inner)]),
      root.edge   = unname(root.edge)
      )

  }

  graph_to_phylo <- function(tree,
                             output_format = output.format,
                             solve_multichotomies = solve.multichotomies){

    g <- tree@tree
    edge_df <- igraph::as_data_frame(g, what = 'edges')
    node_df <- igraph::as_data_frame(g, what = 'vertices')

    edge_matrix <- matrix(cbind(edge_df$from, edge_df$to), ncol = 2)
    weights <- edge_df$weight

    if(any(node_df$node_type == 'germline')){
      if(length(which(node_df$node_type == 'germline')) > 1){
        stop('Cannot create phylo objects for trees with multiple germlines / trees joined with multiple germlines')
      }
      phylo <- to_phylo(x = edge_matrix, root.edge = which(node_df$node_type == 'germline'))
    }else{
      phylo <- to_phylo(x = edge_matrix)
    }

    phylo$edge.length <- weights

    if(solve_multichotomies){
      phylo <- ape::multi2di(phylo)
    }

    if(output_format == 'treedata'){
      requireNamespace('tidytree')
      node_df <- igraph::as_data_frame(g, what = 'vertices')

      tibble_tree <- tidytree::as_tibble(phylo)

      phylo <- merge(tibble_tree, node_df, by = 'label', all.x = TRUE) %>%
               tidytree::as.treedata()
    }
    return(phylo)
  }

  if(inherits(trees, 'list')){
    for(i in 1:length(trees)){
      phylo_list <- vector(mode = "list", length = length(trees[[i]]))

      if(parallel){
        requireNamespace('parallel')
        cores <- parallel::detectCores()
        phylo_list <- parallel::mclapply(trees[[i]], mc.cores = cores,
                                                 FUN = function(x) {x %>% graph_to_phylo()

                                                                    })
      }else{
        phylo_list <- lapply(trees[[i]], function(x) {x %>% graph_to_phylo()
                                                              })
      }
      for(j in 1:length(trees[[i]])){
        trees[[i]][[j]]@phylo <- ape::as.phylo(phylo_list[[j]])
      }

    }

  }else if(inherits(trees, 'AntibodyForests')){
    trees@phylo <- trees %>% graph_to_phylo() %>% ape::as.phylo()


  }else{
    stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
  }

  return(trees)
}
