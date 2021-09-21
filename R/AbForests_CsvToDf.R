#' Convert list of csvs, to nested list of data.frames

#' @description CsvToDf converts a list of csv files, which are clone lineages to a list of data.frames.
#' @param files a list of csv files. Each csv file may contain information concerning full length heavy and light chain sequences, CDRH3 and CDRL3 sequences, the type of isotype and the transcriptional cluster that corresponds to each of these sequences.
#' @return  graphs a list of data.frames. Each data.frame contains 2 columns, one that describes the sequences and the other the type of information (isotype or cluster) is considered in the analysis. All these cases are determined by the user.
#' @export
#' @seealso AntibodyForest
#' @examples
#' \dontrun{
#' CsvToDf(files)
#'}

AbForests_CsvToDf<-function(files){
  graphs<-c()
  i<-0
  for (file in 1:length(files)){
    i<-i+1
    names_col<-names(Filter(is.factor, files[file]))
    data <- data.frame(lapply(files[file], as.character), stringsAsFactors=FALSE)
    temp <- strsplit(names_col, "\\.")
    temp<-gsub('[\"]', '', temp)
    temp<-gsub("^c\\(|\\)$", "", temp)
    names(data)<-temp
    data<-data.frame(data,do.call(rbind,stringr::str_split(data[,1],",")))[,-1]
    temp <- unlist(strsplit(temp, ","))
    colnames(data)<-temp
    data$` Cluster`<-paste0("#",data$` Cluster`)
    row_germ<-unname(which(sapply(data$` Cluster`, function(x) any(grepl("germline", x)))))
    data$` Cluster`[row_germ]<-substring(data$` Cluster`[row_germ], 2)
    data$` Cluster`[row_germ]<-sub("(germline)", "\\1#\\2", data$` Cluster`[row_germ])
    graphs[[i]]<-data

  }
  return(graphs)
}



