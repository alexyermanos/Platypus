#' Split single cell immune repertoire into sub-repertoires by isotype based on number of B cells

#' @description SubRepertoiresByCells separates the single cell immune repertoire into 5 sub-repertoires taking into account the number of cells. The goal is to determine the majority isotype per each network in the immune repertoire. Therefore, each sub-repertoire is dominated by isotype IGG, IGA, IGM, other and if there is an equal number of IGA and IGG isotypes in a network, IGA-IGG category exists respectively. In particular, in case there exists a tie in the number of IGA and IGM, the network is considered to contain IGA as majority isotype, while the same number of IGG and IGM in the network categorize this network as containing IGG as majority isotype. The function receives the output of CsvToDf or original data and can be given as input to ConvertStructure function.
#' @param list a list of data.frames. Each data.frame represents a clone lineage and separates initial input data into subsets of networks.
#' @return  list a nested list of 5 sub-lists of data.frames. Each sub-list corresponds to the set of networks, in which a majority isotype is specifyied. list[[1]] or list$list_IGHG contains the networks, in data.frame format, with more IGG isotypes, list[[2]] or list$list_IGHA contains the networks, in data.frame format, with more IGA isotypes, list[[3]] or list$list_IGHM contains the networks, in data.frame format, with more IGM isotypes, list[[4]] or list$list_IGAG contains the networks, in data.frame format, with a tie in IGA and IGG isotypes and list[[5]] or list$list_other contains the networks, in data.frame format, with other isotypes apart from the aforementioned combinations.
#' @export
#' @seealso ConvertStructure, CsvToDf
#' @examples
#' \dontrun{
#' SubRepertoiresByCells(list)
#'}


AbForests_SubRepertoiresByCells<-function(list){
  if (any(sapply(list,is.list))) {
    cont<-unlist(sapply(list,sapply,function(z) any(grepl(paste(c("VDJ_sequence_nt_trimmed","VJ_sequence_nt_trimmed","VDJ_cgene","VDJ_trimmed_ref","VJ_trimmed_ref","clonotype_id_10x","orig_barcode","seurat_clusters"), collapse = "|"), names(z)))))
    if(any(cont==TRUE)){
      list<-lapply(list,lapply, function(k) .SPLITLINEAGES(k))
      n <- unique(unlist(lapply(list, lapply,names)))
      names(n) <- n
      list<-lapply(n, lapply,function(ni) (lapply(list, lapply,`[[`, ni)))
      list<-lapply(list, lapply,function(x) lapply(x, function(z) Filter(length, z)))
      list<-lapply(list, lapply,function(x) Filter(length, x))
      list<-lapply(list, function(z) unlist(z,recursive=FALSE))
      list<- lapply(list, function(z) unlist(z,recursive=FALSE))
    }else if(any(unlist(sapply(list,function(z) any(grepl(paste(c("X","HCvgene","LCjgene","HCjgene","LCvgene","HCcdr3","LCcdr3"," Cluster","HCisotype","HCLC"), collapse = "|"), names(z))))))==TRUE){
      list<-lapply(list, function(k) .SPLITLINEAGES(k))
      n <- unique(unlist(lapply(list, names)))
      names(n) <- n
      list<-lapply(n, function(ni) (lapply(list, `[[`, ni)))
      list<-lapply(list, function(x) Filter(length, x))

    }else{

      list<-lapply(list,lapply, function(k) .SPLITLINEAGES(k))
      list<-lapply(list, lapply,function(x) utils::head(x))
      list<-lapply(list, function(x) Filter(function(k) length(k)>0, x))
      list <- do.call(Map, c(c,list))
      list<-do.call(Map,c(c,list))
      n <- 2
      list <- lapply(list,function(x) split(x, as.integer(gl(length(x), n,
                                                             length(x))), function(z)
                                                               lapply(z, function(y) data.frame(x=z$Seq,y=z$isotype))))

      list<-lapply(list, function(x) Filter(length, x))
    }
  }
  return(list)

}

