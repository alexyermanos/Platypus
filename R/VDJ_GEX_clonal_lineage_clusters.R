#' only Platypus v2 Integrates the transcriptional cluster information into the clonal lineages. This requires that automate_GEX, VDJ_clonal_lineages, and VDJ_GEX_integrate have already been ran. The transcriptional cluster will be added to the end of the Name for each sequence.
#' @param VDJ_clonal_lineages.output Output from VDJ_clonal_lineages. This should be nested list, with the outer list element corresponding to the individual repertoire and the inner list corresponding to individual clonal lineages based on the initial clonotyping strategy in the form of a dataframe with either Seq or Name. The Name currently contains the barcode following the last "_".
#' @param VDJ_GEX_integrate.output The output from the VDJ_GEX_integrate function that is performed on the VDJ_per_clone level. This involves a nested list where the outer list corresponds to the repertoire and inner lists correspond to specific clones based on the clonotyping strategy.
#' @return a nested list in the identical format to the VDJ_clonal_lineages.output but the name of each sequence will have been changed to include the transcriptional cluster corresponding to that barcode from the GEX library. This requires first running the
#' @export
#' @examples
#' \dontrun{
#' clonal_lineages <- VDJ_clonal_lineages(call_MIXCR.output=call_MIXCR_output
#' , VDJ_extract_germline.output=VDJ_extract_germline_output
#' ,as.nucleotide=FALSE,with.germline=TRUE)
#'}
#'
VDJ_GEX_clonal_lineage_clusters <- function(VDJ_GEX_integrate.output,
                                            VDJ_clonal_lineages.output){
  for(i in 1:length(VDJ_clonal_lineages.output)){
    for(j in 1:length(VDJ_clonal_lineages.output[[i]])){
      barcodes <- gsub(".*_","",VDJ_clonal_lineages.output[[i]][[j]]$Name)
      for(k in 1:length(VDJ_clonal_lineages.output[[i]][[j]]$Name)){
        if(barcodes[k] %in% VDJ_GEX_integrate.output[[i]][[j]]$barcode){
          VDJ_clonal_lineages.output[[i]][[j]]$Name[k] <- paste(VDJ_clonal_lineages.output[[i]][[j]]$Name[k],"_cluster",VDJ_GEX_integrate.output[[i]][[j]]$cluster_membership[which(VDJ_GEX_integrate.output[[i]][[j]]$barcode==barcodes[k])],sep="")
        }
        else{
          VDJ_clonal_lineages.output[[i]][[j]]$Name[k] <- paste(VDJ_clonal_lineages.output[[i]][[j]]$Name[k],"_clusterUnknown",sep="")
        }
      }
    }
  }
  return(VDJ_clonal_lineages.output)
}
