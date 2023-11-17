#'Similarity networks based on CDR3 regions
#'
#' @description Creates a similarity network where clones with similar CDR3s are connected.
#' @param VDJ Either (for platypus version "v2") output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire, OR (for platypus version "v3") the the VDJ matrix output of the VDJ_GEX_matrix() function (VDJ.GEX.matrix.output[[1]])
#' @param distance.cutoff The threshold Levenshtein distance for which two nodes will be connected on the similarity network.
#' @param per.sample logical value indicating if a single networks should be produced for each mouse.
#' @param hcdr3.only logical value indicating if the network is based on heavy chain cdr3s (hcdr3.only = T) or pasted heavy and light chain cdr3s (hcdr3.only = F), works for platypus.version 3 only
#' @param platypus.version Character. Defaults to "v3". Can be "v2" or "v3" dependent on the input format
#' @param known.binders Either a character vector with cdr3s of known binders or a data frame with cdr3s in the first and the corresponding specificity in the second column. If this parameter is defined, the output will be a network with only edges between known binders and the repertoire nodes and edges between the known binders that have at least one edge to a repertoire node
#' @param is.bulk logical value indicating whether the VDJ input was generated from bulk-sequencing data using the bulk_to_vgm function. If is.bulk = T, the VDJ_network function is compatible for use with bulk data. Defaults to False (F).
#' @return returns a list containing networks and network information. If per.sample is set to TRUE then the result will be a network for each repertoire. If per.sample ==F, output[[1]] <- will contain the network, output[[2]] will contain the dataframe with information on each node, such as frequency, mouse origin etc. output[[3]] will contain the connected index - these numbers indicate that the nodes are connected to at least one other node. output[[4]] contains the paired graph - so the graph where only the connected nodes are drawn.
#' @export
#' @examples
#' \dontrun{
#' #Platypus v2
#' #network_out <- VDJ_network(VDJ = VDJ_analyze.out[[1]],per.sample = TRUE,distance.cutoff = 2)
#' #Platypus v3
#' network_out <- VDJ_network(VDJ = Platypus::small_vgm[[1]],per.sample = FALSE,distance.cutoff = 2)
#'}

VDJ_network <- function(VDJ,
                        distance.cutoff,
                        per.sample,
                        platypus.version,
                        known.binders,
                        hcdr3.only,
                        is.bulk){

  connected <- NULL
  Nr_of_VDJ_chains <- NULL
  Nr_of_VJ_chains <- NULL

  if(missing(is.bulk)) is.bulk <- F
  if(missing(platypus.version)) platypus.version <- "v3" ##START V2
  if(missing(distance.cutoff)) distance.cutoff <- 3
  if(missing(known.binders)) known.binders <- F
  if(missing(hcdr3.only)) hcdr3.only <- F

  #Naming compatibility
  VDJ.matrix <- VDJ
  VDJ <- NULL

  #START V2
  if(platypus.version == "v2"){

    clonotype.list <- VDJ.matrix

    if(per.sample==F & known.binders == F){
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
    else if(per.sample==T & known.binders == F){
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
    }  ##END V2
  }
  else if(platypus.version == "v3"){

    if (is.bulk == F) {
      #filtering for max 1VDJ 1VJ chain
      VDJ.matrix <- subset(VDJ.matrix, Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)
    }
    clonotype.list <- VDJ.matrix
    clonotype.list$CDR3_aa_pasted <- paste0(clonotype.list$VDJ_cdr3s_aa, clonotype.list$VJ_cdr3s_aa)

    if(per.sample==F & inherits(known.binders,'logical')){
      if (hcdr3.only == T){
        clonotype.feature <- clonotype.list$VDJ_cdr3s_aa
      } else if (hcdr3.only == F){
        clonotype.feature <- clonotype.list$CDR3_aa_pasted
      }

      distance_matrix <- stringdist::stringdistmatrix(clonotype.feature,clonotype.feature,method = "lv")
      paired_network_aa <- distance_matrix
      clonotype.list$color <- grDevices::rainbow(length(unique(clonotype.list$sample_id)))[clonotype.list$sample_id]
      paired_network_aa[paired_network_aa>distance.cutoff] <- NA
      gc()
      paired_network_aa[paired_network_aa<=distance.cutoff] <- 1
      gc()
      diag(paired_network_aa)<-NA
      gc()
      paired_graph <- igraph::graph_from_adjacency_matrix(paired_network_aa, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
      connected_index <- which(rowSums(paired_network_aa,na.rm = T)>0 & colSums(paired_network_aa,na.rm = T)>0)
      paired_network_aa_connect <- paired_network_aa[connected_index,connected_index]
      out.graph <- igraph::graph_from_adjacency_matrix(paired_network_aa_connect, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
      outlist <- list()
      outlist[[1]] <- out.graph
      outlist[[2]] <- clonotype.list
      outlist[[3]] <- connected_index
      outlist[[4]] <- paired_graph
      return(outlist)
    }
    else if(per.sample==T & inherits(known.binders,'logical')){

      outlist <- list()


      clonotype_list <- list()
      for(i in 1:length(unique(clonotype.list$sample_id))){
        clonotype_list[[i]] <- subset(clonotype.list, clonotype.list$sample_id == unique(clonotype.list$sample_id)[i])
      }
      names(clonotype_list) <- unique(clonotype_list$sample_id)
      message(paste0("Sample order: ", paste0(unique(clonotype.list$sample_id), collapse = " ; ")))

      clonotype.list <- clonotype_list

      graph.list <- list()
      index.list <- list()
      graph.paired.list <- list()

      for(i in 1:length(clonotype.list)){

        if (hcdr3.only == T){
          clonotype.feature <- clonotype.list[[i]]$VDJ_cdr3s_aa
        } else if (hcdr3.only == F){
          clonotype.feature <- clonotype.list[[i]]$CDR3_aa_pasted
        }

        distance_matrix <- stringdist::stringdistmatrix(clonotype.feature,clonotype.feature,method = "lv")
        paired_network_aa <- distance_matrix
        paired_network_aa[paired_network_aa>distance.cutoff] <- NA
        gc()
        paired_network_aa[paired_network_aa<=distance.cutoff] <- 1
        gc()
        diag(paired_network_aa)<-NA
        gc()
        paired_graph <- igraph::graph_from_adjacency_matrix(paired_network_aa, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
        connected_index <- which(rowSums(paired_network_aa,na.rm = T)>0 & colSums(paired_network_aa,na.rm = T)>0)
        paired_network_aa_connect <- paired_network_aa[connected_index,connected_index]
        out.graph <- igraph::graph_from_adjacency_matrix(paired_network_aa_connect, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)

        graph.list[[i]] <- out.graph
        index.list[[i]] <- connected_index
        graph.paired.list[[i]] <- paired_graph
      }
      outlist[[1]] <- graph.list
      outlist[[2]] <- clonotype.list
      outlist[[3]] <- index.list
      outlist[[4]] <- graph.paired.list
      return(outlist)
    }
    else if(per.sample==F & !inherits(known.binders,'logical')){ #STAR GLOBAL known binders

      known.binders <- data.frame(known.binders)
      comb_cdrs <- append(known.binders[,1], clonotype.feature) # append cdrs
      comb_cdrs.df <- data.frame(cdrs = comb_cdrs) # append cdrs
      comb_cdrs.df$isotype <- c(rep("IGHG",nrow(known.binders)),substr(clonotype.list$VDJ_cgene,1,4))
      comb_cdrs.df$isotype[comb_cdrs.df$isotype==''] <- NA
      comb_cdrs.df$isotype[is.na(comb_cdrs.df$isotype)] <- "unknown"
      comb_cdrs.df$mouse <-  c(rep("Known binder",nrow(known.binders)), clonotype.list$sample_id)
      comb_cdrs.df$vgene <-  c(rep("Known binder",nrow(known.binders)), clonotype.list$VDJ_vgene)
      comb_cdrs.df$binder <-  c(rep("Known binder",nrow(known.binders)), rep("Repertoire",nrow(clonotype.list)))
      comb_cdrs.df$frequency <-  c(rep(2, nrow(known.binders)),clonotype.list$new_clonal_frequency)

      if (ncol(known.binders)> 1){
        comb_cdrs.df$specificity = c(known.binders[,2],rep("Repertoire",nrow(clonotype.list)))
      }

      #START matrix calculations
      paired_network_aa <- stringdist::stringdistmatrix(comb_cdrs,comb_cdrs,method = "lv")
      gc()
      paired_network_aa[paired_network_aa>distance.cutoff] <- NA
      gc()
      paired_network_aa[paired_network_aa<=distance.cutoff] <- 1 ##remove edges above distance cutoff
      gc()
      paired_network_aa[(nrow(known.binders)+1):length(comb_cdrs),(nrow(known.binders)+1):length(comb_cdrs)] <-NA #remove edges between repertoire nodes
      gc()
      diag(paired_network_aa)<-NA
      gc()

      delete_index <- which(rowSums(paired_network_aa[1:nrow(known.binders),(nrow(known.binders)+1):length(comb_cdrs)],na.rm = T)==0 &
                              colSums(paired_network_aa[(nrow(known.binders)+1):length(comb_cdrs),1:nrow(known.binders)],na.rm = T)==0) #remove known binders without edges to repertoire nodes
      paired_network_aa[delete_index,1:nrow(known.binders)]<-NA
      paired_network_aa[1:nrow(known.binders),delete_index]<-NA

      connected_index <- which(rowSums(paired_network_aa,na.rm = T)>0 & colSums(paired_network_aa,na.rm = T)>0)
      paired_network_aa_connect <- paired_network_aa[connected_index,connected_index]
      #stop matrix calculations
      out.graph <- igraph::graph_from_adjacency_matrix(paired_network_aa_connect, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
      outlist <- list()
      outlist[[1]] <- out.graph
      outlist[[2]] <- comb_cdrs.df[connected_index,]
      outlist[[3]] <- connected_index
      return(outlist)
    }#STAR GLOBAL known binders

    else if(per.sample==T & !inherits(known.binders,'logical')){ #START per sample known binders

      outlist <- list()


      clonotype_list <- list()
      for(i in 1:length(unique(clonotype.list$sample_id))){
        clonotype_list[[i]] <- subset(clonotype.list, clonotype.list$sample_id == unique(clonotype.list$sample_id)[i])
      }


      names(clonotype_list) <- unique(clonotype_list$sample_id)
      message(paste0("Sample order: ", paste0(unique(clonotype.list$sample_id), collapse = " ; ")))

      clonotype.list <- clonotype_list
      known.binders <- data.frame(known.binders)

      graph.list <- list()
      cdr.list <- list()
      index.list <- list()

      for(i in 1:length(clonotype.list)){
        if (hcdr3.only == T){
          clonotype.feature <- clonotype.list[[i]]$VDJ_cdr3s_aa
        } else if (hcdr3.only == F){
          clonotype.feature <- clonotype.list[[i]]$CDR3_aa_pasted
        }

        comb_cdrs <- append(known.binders[,1], clonotype.feature) # append cdrs
        comb_cdrs.df <- data.frame(cdrs = comb_cdrs) # append cdrs
        comb_cdrs.df$isotype <- c(rep("IGHG",nrow(known.binders)),substr(clonotype.list[[i]]$VDJ_cgene,1,4))
        comb_cdrs.df$isotype[comb_cdrs.df$isotype==''] <- NA
        comb_cdrs.df$isotype[is.na(comb_cdrs.df$isotype)] <- "unknown"
        comb_cdrs.df$mouse <-  c(rep("Known binder",nrow(known.binders)), clonotype.list[[i]]$sample_id)
        comb_cdrs.df$vgene <-  c(rep("Known binder",nrow(known.binders)), clonotype.list[[i]]$VDJ_vgene)
        comb_cdrs.df$binder <-  c(rep("Known binder",nrow(known.binders)), rep("Repertoire",nrow(clonotype.list[[i]])))
        comb_cdrs.df$frequency <-  c(rep(2, nrow(known.binders)),clonotype.list[[i]]$new_clonal_frequency)

        if (ncol(known.binders)> 1){
          comb_cdrs.df$specificity = c(known.binders[,2],rep("Repertoire",nrow(clonotype.list[[i]])))
        }

        #START matrix calculations
        paired_network_aa <- stringdist::stringdistmatrix(comb_cdrs,comb_cdrs,method = "lv")
        gc()
        paired_network_aa[paired_network_aa>distance.cutoff] <- NA
        gc()
        paired_network_aa[paired_network_aa<=distance.cutoff] <- 1 ##remove edges above distance cutoff
        gc()
        paired_network_aa[(nrow(known.binders)+1):length(comb_cdrs),(nrow(known.binders)+1):length(comb_cdrs)] <-NA #remove edges between repertoire nodes
        gc()
        diag(paired_network_aa)<-NA
        gc()

        delete_index <- which(rowSums(paired_network_aa[1:nrow(known.binders),(nrow(known.binders)+1):length(comb_cdrs)],na.rm = T)==0 &
                                colSums(paired_network_aa[(nrow(known.binders)+1):length(comb_cdrs),1:nrow(known.binders)],na.rm = T)==0) #remove known binders without edges to repertoire nodes
        paired_network_aa[delete_index,1:nrow(known.binders)]<-NA
        paired_network_aa[1:nrow(known.binders),delete_index]<-NA

        connected_index <- which(rowSums(paired_network_aa,na.rm = T)>0 & colSums(paired_network_aa,na.rm = T)>0)
        paired_network_aa_connect <- paired_network_aa[connected_index,connected_index]
        #stop matrix calculations
        out.graph <- igraph::graph_from_adjacency_matrix(paired_network_aa_connect, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)

        graph.list[[i]] <- out.graph
        cdr.list[[i]] <- comb_cdrs.df[connected_index,]
        index.list[[i]] <- connected_index
      }#STOP sample loop
    outlist[[1]] <- graph.list
    outlist[[2]] <- cdr.list
    outlist[[3]] <- index.list
    return(outlist)
    }#STOP per sample known binders

  }
}
