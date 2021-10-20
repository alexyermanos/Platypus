#' Platypus V2: Integrates VDJ and gene expression libraries by providing cluster membership seq_per_vdj object and the index of the cell in the Seurat RNA-seq object.
#' @param GEX.object A single seurat object from automate_GEX function. This will likely be supplied as automate_GEX.output[[1]].
#' @param VDJ.per.clone Output from the VDJ_per_clone function. Each element in the list should be found in the output from the automate_GEX function.
#' @return Returns a dataframe containing repertoire information, such as isotype, CDR sequences, mean number of UMIs. This output can be supplied to furhter packages VDJ_extract_sequences and VDJ_GEX_integrate
#' @export
#' @examples
#' \dontrun{
#' GEX_clonotype(GEX.object=automate.GEX.output[[1]], VDJ.per.clone=vdj.per.clone.output)
#'}
#'
GEX_clonotype <- function(GEX.object,VDJ.per.clone){

  GEX.object$clonotype_id <- rep(NA,length(colnames(GEX.object)))
  GEX.object$expanded <- rep(NA,length(colnames(GEX.object)))
  GEX.object$clone_rank <- rep(NA,length(colnames(GEX.object)))

  for(i in 1:length(VDJ.per.clone)){
    holding_bar <- utils::txtProgressBar(min = 0, max = 1, initial = 0, char = "=",
                                  width = NA, style = 1, file = "")
    message(paste(i,"from",length(VDJ.per.clone),"repertoires"))
    for(j in 1:length(VDJ.per.clone[[i]])){
      utils::setTxtProgressBar(value = j/length(VDJ.per.clone[[i]]),pb = holding_bar)
      barcodes <- gsub(VDJ.per.clone[[i]][[j]]$barcode,pattern = "-1",replacement = "")

      for(k in 1:length(barcodes)){
        index <- which(GEX.object$sample_id==i & GEX.object$barcode==barcodes[k])
        GEX.object$clonotype_id[index] <- VDJ.per.clone[[i]][[j]]$clonotype[k]
        GEX.object$expanded[index] <- nrow(VDJ.per.clone[[i]][[j]])
        GEX.object$clone_rank[index] <- j
      }
    }
  }
  return(GEX.object)
}
