#' Creates a similarity network where clones with similar CDR3s are connected.
#' @param clonotype.list dataframe based on the clonotypes.csv file from cellranger vdj output. If a single network is desired then just a subset of the multi-sample clonotype.list should be supplied.
#' @param distance.cutoff The threshold Levenshtein distance for which two nodes will be connected on the similarity network.
#' @param per.mouse logical value indicating if a single networks should be produced for each mouse.
#' @param connected logical value indicating if the connected network should be produced as output. This will result in filtering out all nodes that are not connected in the output network. Only relevant if per.mouse is set to false.
#' @return returns a list containing networks and network information. If per.mouse is set to TRUE then the result will be a network for each repertoire. If per.mouse ==F, output[[1]] <- will contain the network, output[[2]] will contain the dataframe with information on each node, such as frequency, mouse origin etc. output[[3]] will contain the connected index - these numbers indicate that the nodes are connected to at least one other node. output[[4]] contains the paired graph - so the graph where only the connected nodes are drawn.
#' @export
#' @examples
#' \dontrun{
#' VDJ_network(my_VDJ_analyze_final_check[1:1],per.mouse = T,distance.cutoff = 2)
#'}
VDJ_network <- function(clonotype.list,distance.cutoff,per.mouse,connected){
  require(stringdist)
  if(per.mouse==F){
  for(i in 1:length(clonotype.list)) clonotype.list[[i]]$mouse <- rep(i,nrow(clonotype.list[[i]]))
  clonotype_rbind <- do.call("rbind",clonotype.list)
  distance_matrix <- stringdist::stringdistmatrix(clonotype_rbind$CDR3_aa_pasted,clonotype_rbind$CDR3_aa_pasted,method = "lv")
  paired_network_aa <- distance_matrix
  clonotype_rbind$color <- grDevices::rainbow(length(unique(clonotype_rbind$mouse)))[clonotype_rbind$mouse]
  diag(paired_network_aa) <- NA
  paired_network_aa[paired_network_aa<distance.cutoff] <- 1
  paired_network_aa[paired_network_aa>=distance.cutoff] <- 0
  paired_graph <- igraph::graph_from_adjacency_matrix(paired_network_aa, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
  connected_index <- which(rowSums(paired_network_aa,na.rm = T)>0 & colSums(paired_network_aa,na.rm = T)>0)
  paired_network_aa_connect <- paired_network_aa[connected_index,connected_index]
  out.graph <- igraph::graph_from_adjacency_matrix(paired_network_aa_connect, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
  outlist <- list()
  outlist[[1]] <- out.graph
  outlist[[2]] <- clonotype_rbind
  outlist[[3]] <- connected_index
  outlist[[4]] <- paired_graph
  return(outlist)
  }
  else if(per.mouse==T){
    outlist <- list()

    for(i in 1:length(clonotype.list)){
      distance_matrix <- stringdist::stringdistmatrix(clonotype.list[[i]]$CDR3_aa_pasted,clonotype.list[[i]]$CDR3_aa_pasted,method = "lv")
      paired_network_aa <- distance_matrix
      diag(paired_network_aa) <- NA
      paired_network_aa[paired_network_aa<distance.cutoff] <- 1
      paired_network_aa[paired_network_aa>=distance.cutoff] <- 0
      paired_graph <- igraph::graph_from_adjacency_matrix(paired_network_aa, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
      connected_index <- which(rowSums(paired_network_aa,na.rm = T)>0 & colSums(paired_network_aa,na.rm = T)>0)
      paired_network_aa_connect <- paired_network_aa[connected_index,connected_index]
      out.graph <- igraph::graph_from_adjacency_matrix(paired_network_aa_connect, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
      if(connected==TRUE) outlist[[i]] <- out.graph
      else if(connected==FALSE) outlist[[i]] <- paired_graph
    }
    return(outlist)
  }
}
