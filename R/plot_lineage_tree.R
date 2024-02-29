#' Plot B cell lineage tree of clonotype from AntibodyForests object
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description 
#' @param AntibodyForests_object list - AntibodyForests object as obtained from the AntibodyForests() function in Platypus.
#' @param sample string - denotes the sample that contains the clonotype.
#' @param clonotype string - denotes the clonotype from which the B cell lineage tree should be plotted.
#' @return Plots B cell lineage tree for the specified clonotype.
#' @export
#' @examples
#' \dontrun{
#' plot_lineage_tree(AntibodyForests_object,
#'                   sample = "S1",
#'                   clonotype = "clonotype1")
#'}

plot_lineage_tree <- function(AntibodyForests_object,
                              sample,
                              clonotype){
  
  # Retrieve igraph object from AntibodyForests object
  tree <- AntibodyForests_object[[sample]][[clonotype]][[3]]
  
  # Arrange the nodes by using 'germline' node as the root and by directing the tree downwards using the 'igraph::layout_as_tree()' function
  layout <- igraph::layout_as_tree(tree, root = "germline")
  
  # Define node labels
  igraph::V(tree)$label <- ifelse(igraph::V(tree)$name == "germline", "G", 
                                  ifelse(startsWith(igraph::V(tree)$name, "node"), gsub(pattern = "node", replacement = "", igraph::V(tree)$name), igraph::V(tree)$name))
  
  # Define node colors (the germline node is colored orange, the recovered sequence nodes are colored lightblue, and the unrecovered (internal) sequence nodes are colored grey)
  igraph::V(tree)$color <- ifelse(igraph::V(tree)$label == "G", "orange",
                                  ifelse(startsWith(igraph::V(tree)$name, "node"), "lightblue", "grey"))
  
  # Plot tree
  igraph::plot.igraph(tree, layout = layout,
                      vertex.label.dist = 0,
                      vertex.size = 10, 
                      vertex.label.cex = 0.8,
                      edge.arrow.size = 0.1)
}