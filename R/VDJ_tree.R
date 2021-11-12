#' Produces neighbor joining phylogenetic trees from the output of VDJ_clonal_lineages
#' @param clonal.lineages Output from VDJ_clonal_lineages. This should be nested list, with the outer list element corresponding to the individual repertoire and the inner list corresponding to individual clonal lineages based on the initial clonotyping strategy in the form of a dataframe with either Seq or Name.
#' @param molecular.marker Biological observation on which phylogenetic trees will be based on. Default is set to hc.lc.full.length. Options: hc.lc.full.length, hc.lc.full.length, hc.cdr3s,hc.lc.cdr3s
#' @param with.germline Logical specifying if the germline should be set as outgroup (germline appears as an additional node). Default is set to TRUE.
#' @param min.sequences integer value specifying the minimum number of sequences to be allowed for clonal lineages. Default is 3.
#' @param max.sequences integer value specifying the maximum number of sequences to be allowed for clonal lineages. Default is 500
#' @param normalize.germline.length Logical determining whether or not the branch length separating the germline from the first internal node should be normalized. Potentially useful for visualization if the remainder tips are far from the root. Default is TRUE.
#' @param unique.sequences Logical indicating if those cells containing identical VDJRegion sequences should be merged into single nodes and have their variant added as the tip label. Default is TRUE.
#' @param platypus.version  Only platypus "v3" is supported. 
#' @return Returns a nested list of phylogenetic trees. The output[[i]][[j]] corresponds to the j'th clone in the i'th input repertoire. plot(output[[i]][[j]]) should display the phylogenetic tree if the ape package is loaded.
#' @export
#' @examples
#' \dontrun{
#' vdj.tree <- VDJ_tree(clonal.lineages = VDJ.clonal.lineage.output,with.germline=T,min.sequences = 5,max.sequences = 30,unique.sequences = T)
#'}
VDJ_tree <- function(clonal.lineages,molecular.marker, with.germline,min.sequences,max.sequences,normalize.germline.length,unique.sequences,platypus.version){
  require(ape)
  require(phytools)
  require(stringdist)
  if(missing(max.sequences)) max.sequences <- 500
  if(missing(min.sequences)) min.sequences <- 3
  if(missing(unique.sequences)) unique.sequences <- T
  if(missing(with.germline)) with.germline <- T
  if(missing(normalize.germline.length)) normalize.germline.length <- T
  if(missing(molecular.marker)) molecular.marker <- "hc.lc.full.length"
  if(missing(platypus.version)) platypus.version <- "v3"
  output.tree <- list()
  
  if (platypus.version == "v2"){ ##START v2
    for(i in 1:length(clonal.lineages)){
      output.tree[[i]] <- list()
      for(j in 1:length(clonal.lineages[[i]])){
        if(nrow(clonal.lineages[[i]][[j]])>=min.sequences & nrow(clonal.lineages[[i]][[j]]<=max.sequences)){
          tempdistance <- stringdist::stringdistmatrix(clonal.lineages[[i]][[j]]$Seq,clonal.lineages[[i]][[j]]$Seq,method="lv")
          output.tree[[i]][[j]] <- ape::nj(tempdistance)
          if(with.germline==T){
            output.tree[[i]][[j]] <- phytools::reroot(output.tree[[i]][[j]],node.number = which(grepl(clonal.lineages[[i]][[j]]$Name,pattern = "germline")))
            if(normalize.germline.length==T) output.tree[[i]][[j]]$edge.length[1] <- .5
          }
        }
      }
      print(i)
    }
    if(unique.sequences==T){
      for(i in 1:length(clonal.lineages)){
        for(j in 1:length(clonal.lineages[[i]])){
          dupicated_which <- which(duplicated(clonal.lineages[[i]][[j]]$Seq)==F)
          clonal.lineages.unique <- clonal.lineages[[i]][[j]][dupicated_which,]
          if(nrow(clonal.lineages.unique)>=min.sequences & nrow(clonal.lineages.unique<=max.sequences)){
            clonal.lineages.unique$freq_name <- rep("",nrow(clonal.lineages.unique))
            for(k in 1:nrow(clonal.lineages.unique)){
              clonal.lineages.unique$freq_name[k] <- paste("Seq_",k,"_Freq_",length(which(clonal.lineages[[i]][[j]]$Seq==clonal.lineages.unique$Seq[k])),sep="")
            }
            clonal.lineages.unique.dist <- stringdist::stringdistmatrix(clonal.lineages.unique$Seq,clonal.lineages.unique$Seq)
            rownames(clonal.lineages.unique.dist) <- clonal.lineages.unique$freq_name
            colnames(clonal.lineages.unique.dist) <- clonal.lineages.unique$freq_name
            output.tree[[i]][[j]] <- ape::nj(clonal.lineages.unique.dist)
            output.tree[[i]][[j]] <- phytools::reroot(output.tree[[i]][[j]],node.number = which(grepl(clonal.lineages.unique$Name,pattern = "germline")))
            output.tree[[i]][[j]]$tip.label[which(grepl(clonal.lineages.unique$Name,pattern = "germline"))]
          }
        }
      }
    }
    return(output.tree)
  }##STOP v2
  
  if (platypus.version == "v3"){ ##START v3
  
    if (molecular.marker != "hc.full.length" & molecular.marker != "hc.lc.full.length" & with.germline == T){
      with.germline <- F
      print("The parameter with.germline will be set to false (germline information only available if hc.full.length or hc.lc.full.length selected as molecular.marker)")
    }
    
    if (min.sequences < 3){
      min.sequences <- 3
      print("At least 3 observations needed to build unrooted tree (min.sequences will be set to 3)")
    }
    if (max.sequences >500){
      max.sequences <- 500
      print("At max.500 observations can be used to build unrooted tree (max.sequences will be set to 500)")
    }
    

    #slected molecular marker is saved under clonal.lineages[[i]][[j]]$seq
    if (molecular.marker == "hc.full.length"){  
      for (i in 1:length(clonal.lineages)){
        for (j in 1:length(clonal.lineages[[i]])){
          clonal.lineages[[i]][[j]]$seq <- clonal.lineages[[i]][[j]]$VDJ_sequence_nt_trimmed
        }
      }
    } else if (molecular.marker == "hc.lc.full.length"){  
      for (i in 1:length(clonal.lineages)){
        for (j in 1:length(clonal.lineages[[i]])){
          clonal.lineages[[i]][[j]]$seq <- paste0(clonal.lineages[[i]][[j]]$VDJ_sequence_nt_trimmed,clonal.lineages[[i]][[j]]$VJ_sequence_nt_trimmed)
        }
      }  
    } else if (molecular.marker == "hc.cdr3s"){  
      for (i in 1:length(clonal.lineages)){
        for (j in 1:length(clonal.lineages[[i]])){
          clonal.lineages[[i]][[j]]$seq <- clonal.lineages[[i]][[j]]$VDJ_cdr3s_nt
        }
      }
    } else if (molecular.marker == "hc.lc.cdr3s"){  
      for (i in 1:length(clonal.lineages)){
        for (j in 1:length(clonal.lineages[[i]])){
          clonal.lineages[[i]][[j]]$seq <- paste0(clonal.lineages[[i]][[j]]$VDJ_cdr3s_nt,clonal.lineages[[i]][[j]]$VJ_cdr3s_nt)
        }
      }
    }##stop selecting molecular.marker
    
    if(unique.sequences==F){
      for(i in 1:length(clonal.lineages)){
        output.tree[[i]] <- list()
        for(j in 1:length(clonal.lineages[[i]])){
          if (molecular.marker == "hc.full.length" & with.germline == T){
            germline.node <- names(sort(table(clonal.lineages[[i]][[j]]$VDJ_trimmed_ref),decreasing = T))[1]
          } else if (molecular.marker == "hc.lc.full.length" & with.germline == T){
            germline.node <- names(sort(table(paste0(clonal.lineages[[i]][[j]]$VDJ_trimmed_ref,clonal.lineages[[i]][[j]]$VJ_trimmed_ref)),decreasing = T))[1]
          }
          if(nrow(clonal.lineages[[i]][[j]])>=min.sequences & nrow(clonal.lineages[[i]][[j]])<=max.sequences){ #check if minimal number of cells needed for tree construction present
            clonal.lineages[[i]][[j]]$freq_name <- rep("",nrow(clonal.lineages[[i]][[j]]))
            for(k in 1:nrow(clonal.lineages[[i]][[j]])){
              clonal.lineages[[i]][[j]]$freq_name[k] <- paste("Seq_",k,sep="")
            }
            if (with.germline == T){
              tempdistance <- stringdist::stringdistmatrix(c(clonal.lineages[[i]][[j]]$seq,germline.node),c(clonal.lineages[[i]][[j]]$seq,germline.node))
              rownames(tempdistance) <- c(clonal.lineages[[i]][[j]]$freq_name,"germline")
              colnames(tempdistance) <- c(clonal.lineages[[i]][[j]]$freq_name,"germline")
              output.tree[[i]][[j]] <- ape::nj(tempdistance)
              output.tree[[i]][[j]]<- phytools::reroot(output.tree[[i]][[j]],node.number = length(output.tree[[i]][[j]]$tip.label))
              if(normalize.germline.length==T) output.tree[[i]][[j]]$edge.length[1] <- .5
            }else if (with.germline == F){
              tempdistance <- stringdist::stringdistmatrix(clonal.lineages[[i]][[j]]$seq,clonal.lineages[[i]][[j]]$seq)
              rownames(tempdistance) <- clonal.lineages[[i]][[j]]$freq_name
              colnames(tempdistance) <- clonal.lineages[[i]][[j]]$freq_name
              output.tree[[i]][[j]] <- ape::nj(tempdistance)
            }
          }else {
            output.tree[[i]][[j]] <- paste0("# unique observations: ",length(unique(clonal.lineages[[i]][[j]]$seq)))
          }
        }##stop iterating over clonotypes
      }##stop iterating over repertoires
    }## Stop unique.sequences==F loop
    if(unique.sequences==T){
      for(i in 1:length(clonal.lineages)){
        output.tree[[i]] <- list()
        for(j in 1:length(clonal.lineages[[i]])){
          if (molecular.marker == "hc.full.length" & with.germline == T){
            germline.node <- names(sort(table(clonal.lineages[[i]][[j]]$VDJ_trimmed_ref),decreasing = T))[1]
          } else if (molecular.marker == "hc.lc.full.length" & with.germline == T){
            germline.node <- names(sort(table(paste0(clonal.lineages[[i]][[j]]$VDJ_trimmed_ref,clonal.lineages[[i]][[j]]$VJ_trimmed_ref)),decreasing = T))[1]
          }
          dupicated_which <- which(duplicated(clonal.lineages[[i]][[j]]$seq)==F)
          clonal.lineages.unique <- clonal.lineages[[i]][[j]][dupicated_which,]
          if(nrow(clonal.lineages.unique)>=min.sequences & nrow(clonal.lineages.unique)<=max.sequences){ #check if minimal number of cells needed for tree construction present
            clonal.lineages.unique$freq_name <- rep("",nrow(clonal.lineages.unique))
            for(k in 1:nrow(clonal.lineages.unique)){
              clonal.lineages.unique$freq_name[k] <- paste("Seq_",k,"_Freq_",length(which(clonal.lineages[[i]][[j]]$seq==clonal.lineages.unique$seq[k])),sep="")
            }
            if (with.germline == T){
              clonal.lineages.unique.dist <- stringdist::stringdistmatrix(c(clonal.lineages.unique$seq,germline.node),c(clonal.lineages.unique$seq,germline.node))
              rownames(clonal.lineages.unique.dist) <- c(clonal.lineages.unique$freq_name,"germline")
              colnames(clonal.lineages.unique.dist) <- c(clonal.lineages.unique$freq_name,"germline")
              output.tree[[i]][[j]] <- ape::nj(clonal.lineages.unique.dist)
              output.tree[[i]][[j]]<- phytools::reroot(output.tree[[i]][[j]],node.number = length(output.tree[[i]][[j]]$tip.label))
              if(normalize.germline.length==T) output.tree[[i]][[j]]$edge.length[1] <- .5
            }else if (with.germline == F){
              clonal.lineages.unique.dist <- stringdist::stringdistmatrix(clonal.lineages.unique$seq,clonal.lineages.unique$seq)
              rownames(clonal.lineages.unique.dist) <- clonal.lineages.unique$freq_name
              colnames(clonal.lineages.unique.dist) <- clonal.lineages.unique$freq_name
              output.tree[[i]][[j]] <- ape::nj(clonal.lineages.unique.dist)
            }
          }else {
            output.tree[[i]][[j]] <- paste0("# unique observations: ",length(unique(clonal.lineages.unique$seq)))
          }
        }##stop iterating over clonotypes
      }##stop iterating over repertoires
    }##stop unique.sequences==T loop
    return(output.tree)
  } #STOP v3
}  










