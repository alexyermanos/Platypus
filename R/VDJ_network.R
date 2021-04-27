#' Creates a similarity network where clones with similar CDR3s are connected.
#' @param VDJ.matrix Either (for platypus version "v2") output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire, OR (for platypus version "v3") the the VDJ matrix output of the VDJ_GEX_matrix() function (normally VDJ.GEX.matrix.output[[1]]) 
#' @param distance.cutoff The threshold Levenshtein distance for which two nodes will be connected on the similarity network.
#' @param per.sample logical value indicating if a single networks should be produced for each mouse.
#' @param connected logical value indicating if the connected network should be produced as output. This will result in filtering out all nodes that are not connected in the output network. Only relevant if per.sample is set to false.
#' @param platypus.version Character. Defaults to "v2". Can be "v2" or "v3" dependent on the input format 
#' @return returns a list containing networks and network information. If per.sample is set to TRUE then the result will be a network for each repertoire. If per.sample ==F, output[[1]] <- will contain the network, output[[2]] will contain the dataframe with information on each node, such as frequency, mouse origin etc. output[[3]] will contain the connected index - these numbers indicate that the nodes are connected to at least one other node. output[[4]] contains the paired graph - so the graph where only the connected nodes are drawn.
#' @export
#' @examples
#' \dontrun{
#' VDJ_network(my_VDJ_analyze_final_check[1:1],per.sample = T,distance.cutoff = 2)
#'}
VDJ_network <- function(VDJ.matrix,distance.cutoff,per.sample,connected,platypus.version){
  require(stringdist)
  
  if(missing(platypus.version)) platypus.version <- "v2"
  
  if(platypus.version == "v2"){
  
  clonotype.list <- VDJ.matrix
    
  if(per.sample==F){
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
  else if(per.sample==T){
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
  } else if(platypus.version == "v3"){
    
    #filtering for max 1VDJ 1VJ chain
    VDJ.matrix <- subset(VDJ.matrix, Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)
    
    clonotype.list <- VDJ.matrix
    
    clonotype.list$CDR3_aa_pasted <- paste0(clonotype.list$VDJ_cdr3s_aa, clonotype.list$VJ_cdr3s_aa)
    clonotype.list$mouse <- clonotype.list$sample_id
    
    if(per.sample==F){
      #for(i in 1:length(clonotype.list)) clonotype.list[[i]]$mouse <- rep(i,nrow(clonotype.list[[i]]))
      #clonotype_rbind <- do.call("rbind",clonotype.list)
      clonotype_rbind <- clonotype.list
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
    else if(per.sample==T){
      
      outlist <- list()
      

      clonotype_list <- list()
      for(i in 1:length(unique(clonotype.list$sample_id))){
        clonotype_list[[i]] <- subset(clonotype.list, clonotype.list$sample_id == unique(clonotype.list$sample_id)[i])
      }
      names(clonotype_list) <- unique(clonotype_list$sample_id)
      print(paste0("Sample order: ", paste0(unique(clonotype.list$sample_id), collapse = " ; ")))
      
      clonotype.list <- clonotype_list
      
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
}
