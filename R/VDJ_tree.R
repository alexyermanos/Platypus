#' Produces neighbor joining phylogenetic trees from the output of VDJ_clonal_lineages
#' @param clonal.lineages Output from VDJ_clonal_lineages. This should be nested list, with the outer list element corresponding to the individual repertoire and the inner list corresponding to individual clonal lineages based on the initial clonotyping strategy in the form of a dataframe with either Seq or Name.
#' @param with.germline Logical specifying if the germline should be set as outgroup. Default is set to TRUE.
#' @param min.sequences integer value specifying the minimum number of sequences to be allowed for clonal lineages. Default is 3.
#' @param max.sequences integer value specifying the maximum number of sequences to be allowed for clonal lineages. Default is 500
#' @param normalize.germline.length Logical determining whether or not the branch length separating the germline from the first internal node should be normalized. Potentially useful for visualization if the remainder tips are far from the root. Default is TRUE.
#' @param unique.sequences Logical indicating if those cells containing identical VDJRegion sequences should be merged into single nodes and have their variant added as the tip label. Default is TRUE.
#' @return Returns a nested list of phylogenetic trees. The output[[i]][[j]] corresponds to the j'th clone in the i'th input repertoire. plot(output[[i]][[j]]) should display the phylogenetic tree if the ape package is loaded.
#' @export
#' @examples
#' \dontrun{
#' vdj.tree <- VDJ_tree(clonal.lineages = VDJ.clonal.lineage.output,with.germline=T,min.sequences = 5,max.sequences = 30,unique.sequences = T)
#'}
VDJ_tree <- function(clonal.lineages,with.germline,min.sequences,max.sequences,normalize.germline.length,unique.sequences){
  require(ape)
  require(phytools)
  require(stringdist)
  if(missing(max.sequences)) max.sequences <- 500
  if(missing(max.sequences)) min.sequences <- 3
  if(missing(unique.sequences)) unique.sequences <- T
  if(missing(with.germline)) with.germline <- T
  if(missing(normalize.germline.length)) normalize.germline.length <- T
  output.tree <- list()
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
}
