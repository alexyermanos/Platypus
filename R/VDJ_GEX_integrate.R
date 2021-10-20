#' only Platypus v2 Integrates VDJ and gene expression libraries by providing cluster membership seq_per_vdj object and the index of the cell in the Seurat RNA-seq object.
#' @param GEX.object A single seurat object from automate_GEX function. This will likely be supplied as automate_GEX.output[[1]].
#' @param clonotype.list  Output from either VDJ_analyze or VDJ_clonotype functions. This list should correspond to a single GEX.list object, in which each list element in clonotype.list is found in the GEX.object. Furthermore, these repertoires should be found in the automate_GEX library.
#' @param VDJ.per.clone Output from the VDJ_per_clone function. Each element in the list should be found in the output from the automate_GEX function.
#' @param clonotype.level Logical specifying whether the integration should occur on the cellular level (VDJ_per_clone) or on the clonotype level (e.g. output from VDJ_analyze or VDJ_clonotype). TRUE specifies that the clonotype level will be selected - e.g. the clonotype.list object will now contain information from the GEX object regarding clonal membership.
#' @return Returns a nested list containing information corresponding to either the clonal level or the sequence level, depending on the input argument "clonotype.level". This function essentially will update the output of the analyze_VDJ or the VDJ_per_clone functions.
#' @export
#' @examples
#' \dontrun{
#' testing_integrate <- VDJ_GEX_integrate(GEX.object = automate.gex.output[[1]]
#' ,clonotype.list =  VDJ.analyze.output
#' ,VDJ.per.clone = VDJ.per.clone.output,clonotype.level = TRUE)
#'}

VDJ_GEX_integrate <- function(GEX.object, clonotype.list,VDJ.per.clone,clonotype.level){
  if(clonotype.level==TRUE){
    if(!missing(VDJ.per.clone)) VDJ.per.clone <- NULL
    for(i in 1:length(clonotype.list)){
      clonotype.list[[i]]$majority_cluster <- rep("",nrow(clonotype.list[[i]]))
      clonotype.list[[i]]$cluster_membership_percent <- rep("",nrow(clonotype.list[[i]]))
      clonotype.list[[i]]$cell_index <- rep("",nrow(clonotype.list[[i]]))

      for(j in 1:nrow(clonotype.list[[i]])){

        barcodes.split <- stringr::str_split(string = clonotype.list[[i]]$barcodes[j],pattern = ";")[[1]]
        barcodes.split <- gsub(barcodes.split, pattern = "-1",replacement = "")
        tempindex <- list()
        for(k in 1:length(barcodes.split)){
          tempindex[[k]] <- which(GEX.object$barcode==barcodes.split[k])
        }
        tempindex <- unlist(tempindex)
        clonotype.list[[i]]$cell_index[j] <- gsub(toString(tempindex),pattern = ", ",replacement = ";")
        clonotype.list[[i]]$cluster_membership_percent[j] <- toString(100*table(GEX.object$seurat_clusters[tempindex])/sum(table(GEX.object$seurat_clusters[tempindex])))
        clonotype.list[[i]]$majority_cluster[j] <- names(which.max(table(GEX.object$seurat_clusters[tempindex])))
      }
    }
  }
  else if(clonotype.level==FALSE){
    if(!missing(clonotype.list)) clonotype.list <- NULL
    for(i in 1:length(VDJ.per.clone)){
      for(j in 1:length(VDJ.per.clone[[i]])){
        VDJ.per.clone[[i]][[j]]$cluster_membership <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$cell_index <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        for(k in 1:nrow(VDJ.per.clone[[i]][[j]])){
          if(gsub(VDJ.per.clone[[i]][[j]]$barcode[k],pattern = "-1",replacement = "") %in% colnames(GEX.object)){
            tempindex <- which(GEX.object$barcode==gsub(VDJ.per.clone[[i]][[j]]$barcode[k],pattern = "-1",replacement = "") & GEX.object$sample_id==i)[1]
            VDJ.per.clone[[i]][[j]]$cluster_membership[k] <- GEX.object$seurat_clusters[tempindex]
            VDJ.per.clone[[i]][[j]]$cell_index[k] <- tempindex
          }

        }
      }
    }
  }
  if(clonotype.level==TRUE) return(clonotype.list)
  else if(clonotype.level==F) return(VDJ.per.clone)
}
