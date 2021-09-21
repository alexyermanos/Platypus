#' Conversion of phylogenetic tree to distance matrix

#' @description PlyloToMatrix converts a previously existing phylogenetic tree to a corresponding distance matrix using the cophenetic distance. Then, there is the option to utilize this custom distance matrix as an input distance matrix to AntibodyForest function. The user is responsible for specifying a correct and valid distance matrix. In particular, the size of distance matrix must match the number of sequences for each network in each repertoire.
#' @param tree_name  a plylogenetic tree (phylo object).
#' @return dist_mat The corresponding distance matrix uses the cophenetic distance between two observations that have been clustered. This distance is defined to be the intergroup dissimilarity at which the two observations are first combined into a single cluster.
#' @export
#' @seealso AntibodyForest
#' @examples
#' \dontrun{
#' PlyloToMatrix(tree_name)
#'}
AbForests_PlyloToMatrix<-function(tree_name){
  tree<-ape::read.tree(tree_name)
  dist_mat<-stats::cophenetic(tree)
  return(dist_mat)
}
