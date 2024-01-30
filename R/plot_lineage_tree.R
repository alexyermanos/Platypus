#' Plot B cell lineage tree from AntibodyForests object
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description 
#' @param AntibodyForests_object list - AntibodyForests object as obtained from the AntibodyForests() function in Platypus.
#' @param sample string - denotes the sample that contains the clonotype.
#' @param clonotype string - denotes the clonotype from which the B cell lineage tree should be plotted.
#' @return Plots B cell lineage tree for the specified clonotype.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_plot(AntibodyForests_object,
#'                      sample = "S1",
#'                      clonotype = "clonotype4")
#'}

plot_lineage_tree <- function(AntibodyForests_object,
                              sample,
                              clone){
  
  # Retrieve igraph object from AntibodyForests object
  tree <- AntibodyForests_object[[sample]][[clone]][["lineage.tree"]]
  
  # Arrange the nodes by using 'germline' node as the root and by directing the tree downwards
  layout <- igraph::layout_as_tree(tree, root = which(igraph::V(tree)$node_type == "germline"), circular = FALSE)
  
  # Plot tree
  igraph::plot.igraph(tree, layout = layout, 
                      vertex.label.dist = 0, 
                      edge.arrow.size = 0.1, 
                      vertex.size = 10, 
                      vertex.label.cex = 0.5)
}