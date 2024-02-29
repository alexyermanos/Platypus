#'Converts an igraph network into a phylogenetic tree as a phylo object.

#'@description Converts an igraph network into a phylogenetic tree as a phylo object.

#' @param tree igraph object
#' @param solve.multichotomies boolean - whether to remove multichotomies in the resulting phylogenetic tree using ape::multi2di
#' @return phylogenetic tree
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_phylo(trees)
#'}


AntibodyForests_phylo <- function(tree,
                           solve_multichotomies){
  
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
  
  #Transforms an edge matrix into an object of class phylo
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
    
    #Create phylo object (derived from aphylo package)
    phylo_object(
      edge        = unname(x),
      edge.length = unname(edge.length),
      tip.label   = unname(label[leaves]),
      Nnode       = length(inner) + 1L,
      node.label  = unname(label[c(roots, inner)]),
      root.edge   = unname(root.edge)
    )
    
  }
  
  #Create edge matrix of the igraph tree
  g <- tree
  edge_df <- igraph::as_data_frame(g, what = 'edges')
  node_df <- igraph::as_data_frame(g, what = 'vertices')
  edge_matrix <- matrix(cbind(edge_df$from, edge_df$to), ncol = 2)
  
  #Transform the edge matrix into a phylo tree
  #If there is a node called "germline", root the tree on the germline, otherwise don't root the tree
  if(any(node_df$name == 'germline')){
    if(length(which(node_df$name == 'germline')) > 1){
      stop('Cannot create phylo objects for trees with multiple germlines / trees joined with multiple germlines')
    }
    phylo <- to_phylo(x = edge_matrix, root.edge = which(node_df$name == 'germline'))
  }else{
    phylo <- to_phylo(x = edge_matrix)
  }
  
  #Add edge lengths as branch length to the phylo object
  phylo$edge.length <- edge_df$edge.length
  
  #Transform polytomies into series of bifurcations
  if(solve_multichotomies){
    phylo <- ape::multi2di(phylo)
  }
  
  #Return phylogenetic tree
  return(phylo)
}