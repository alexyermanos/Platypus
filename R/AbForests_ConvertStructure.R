#' Extract transcriptome/isotype information and B cell receptor sequences from single cell immune repertoire formatted as list of data.frames

#' @description ConvertStructure alters a list of clone lineages, represented as data.frames and recovers the type of isotypes or transcriptional clusters and antibody sequences from these clone lineages. It can receive as input the original data or the output of SubRepertoiresByUniqueSeq or SubRepertoiresByCells. Then, the output list is used as input to the RemoveNets or AntibodyForest function.
#' @param list a list of data.frames. Each data.frame may contain information concerning full length heavy and light chain sequences, CDRH3 and CDRL3 sequences, the type of isotype and the transcriptional cluster that corresponds to each of these sequences.
#' @param opt a string with options "isotype" and "cluster". The option "isotype" is utilized when the user desires to do an isotype analysis, while the selection of "cluster" denotes that an analysis based on transcriptome is requested.
#' @param cdr3 variable with values 0 if the user desires to select full length sequences (only when the input is a list of csv files), 1 for sequences in the CDR3 only (only when the input is a list of csv files) and NULL otherwise.
#' @return  list a list of data.frames. Each data.frame contains 2 columns, one that describes the sequences and the other which type of information (isotype or cluster) is considered in the analysis. All these cases are determined by the user.
#' @export
#' @seealso RemoveNets, AntibodyForest
#' @examples
#' \dontrun{
#' ConvertStructure(list,opt="cluster",cdr3=NULL)
#' ConvertStructure(list,opt="isotype",1)
#'}

AbForests_ConvertStructure<-function(list,opt,cdr3){
  if (any(sapply(list, class) == "list")) {
    if(length(cdr3)>0){
      list<-lapply(list, lapply, function(k) .MODIFY_INPUT(k,.REFERENCE_INFO()$type_iso,cdr3,opt))
    }else{
      list<-lapply(list, lapply, function(k) .MODIFY_INPUT(k,.REFERENCE_INFO()$type_iso,NULL,opt))
    }
  }else{
    final_list<-c()
    if(length(cdr3)>0){
      list<-lapply(list, function(k) .MODIFY_INPUT(k,.REFERENCE_INFO()$type_iso,cdr3,opt))
    }else{
      list<-lapply(list, function(k) .MODIFY_INPUT(k,.REFERENCE_INFO()$type_iso,NULL,opt))
    }
    final_list[[1]]<-list
    list<-final_list
  }
  list<-lapply(list, function(x) Filter(function(k) length(k)>0, x))
  return(list)
}


