#' Organizes and extracts full-length sequences for clonal lineage inference. The output sequence can either contain the germline sequence as determined by cellranger or can just contain the sequences contained in each clonal family.
#' @param call_MIXCR.output The output from the call_MIXCR function containing the full-length VDJRegion sequences.
#' @param VDJ_extract_germline.output The output from the VDJ_extract_germline function. This should have the germline information. This needs to be supplied if the with.germline argument is set to true.
#' @param as.nucleotide Logical determining whether the full-length VDJRegion sequence should use nucleotide seqeunce. TRUE indicates nucleotide sequences and FALSE will extract amino acid sequences.
#' @param with.germline Logical determining whether the germline sequence as determined by cellranger should be included in the output list of sequences. If so, the germline will be added to the last row of each dataframe object.
#' @return returns a list containing the sequences for each clonal family as determined by the input clonotyping strategy to call_MIXCR and VDJ_extract_germline. The outer list corresponds to distinct repertoires supplied to the call_MIXCR function (e.g. VDJ.clonal.lineage.output[[i]][[j]] will contain a dataframe of the j'th clone in the i'th repertoire)
#' @export
#' @examples
#' \dontrun{
#' clonal_lineages <- VDJ_clonal_lineages(call_MIXCR.output=call_MIXCR_output, VDJ_extract_germline.output=VDJ_extract_germline_output,as.nucleotide=F,with.germline=T)
#'}
#'
VDJ_clonal_lineages <- function(call_MIXCR.output,VDJ_extract_germline.output,as.nucleotide,with.germline){
  if(missing(VDJ_extract_germline.output)) VDJ_extract_germline.output <- list()
  clonal.lineage.output <- list()
  for(i in 1:length(call_MIXCR.output)){
    clonal.lineage.output[[i]] <- list()
    for(j in 1:length(call_MIXCR.output[[i]])){
      clonal.lineage.output[[i]][[j]] <- data.frame(Seq=rep("",nrow(call_MIXCR.output[[i]][[j]])),
                                                    Name=rep("",nrow(call_MIXCR.output[[i]][[j]])))
      clonal.lineage.output[[i]][[j]]$Name <- paste(call_MIXCR.output[[i]][[j]]$clonotype_id,j,call_MIXCR.output[[i]][[j]]$isotype_hc,call_MIXCR.output[[i]][[j]]$barcode,sep="_")

      if(as.nucleotide==T){
        clonal.lineage.output[[i]][[j]]$Seq <- gsub(call_MIXCR.output[[i]][[j]]$VDJ.NT.HC.LC,pattern = "_",replacement = "")
      }
      else if(as.nucleotide==F){
        clonal.lineage.output[[i]][[j]]$Seq <- gsub(call_MIXCR.output[[i]][[j]]$VDJ.AA.HC.LC,pattern = "_",replacement = "")
      }
      if(with.germline==T){
        if(nrow(VDJ_extract_germline.output[[i]][[1]])==length(call_MIXCR.output[[i]])){
          if(as.nucleotide==T){
            temp_df <- data.frame(Seq=VDJ_extract_germline.output[[i]][[1]]$VDJ.NT.HC.LC[j],Name="germline")
          }
          else if (as.nucleotide==F){
            temp_df <- data.frame(Seq=VDJ_extract_germline.output[[i]][[1]]$VDJ.AA.HC.LC[j],Name="germline")

          }
          clonal.lineage.output[[i]][[j]] <- rbind(clonal.lineage.output[[i]][[j]],temp_df)
        }
      }
      }
  }
  return(clonal.lineage.output)
}

