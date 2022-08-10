

#' This PlatypusML_feature_extraction function takes as input specified features from the first output of the VDJ_GEX_matrix function and encodes 
#' according to the specified strategy. The function returns a matrix containing the encoded extracted features in the order specified in 
#' the input as columns and the different cells as rows. This function should be called as a first step in the process of modeling the VGM data
#' using machine learning. 
#' @param VGM VGM output of the VDJ_GEX_matrix function 
#' @param which.features String vector. Information on which columns of the VDJ input the function should encode
#' @param which.encoding String vector of size 2. Defaults to 'onehot'. Information on which encoding strategy to be used for the two types of sequences:
#' the first entry of the vector corresponds to the nucleotide type of encoding and the second one to the amino acid type of encoding. 
#' If one type of sequence is not among the  
#' Other possible values for amino acid sequences are 'kmer', 'blosum', 'dc', 'tc' or 'topoPCA' and for nucleotide sequences 'kmer'.
#' @param encoding.level String. Specifies on which level the features will be extracted. There are three possible options: "cell" (all available),
#' "clone" (one unique sample per clone), "unique.sequence" (selecting only unique sequences based on a specified sequence (int he unique.sequence argument)). 
#' It defaults to cell.
#' @param unique.sequence String. Needs to be specified only when encoding.level is set to "unique.sequence". The name of the sequence on which unique selection should be 
#' based on.
#' @param parameters.encoding.nt List. Parameters to be used for encoding, if the chosen encoding requires it.
#' 'onehot' -> no parameters necessary, defaults to NULL
#' 'kmer' -> one parameter necessary to set the length of the subsequence, defaults to 3
#' @param parameters.encoding.aa List. Parameters to be used for encoding, if the chosen encoding requires it.
#' 'onehot', 'dc', 'tc' -> no parameters necessary, defaults to NULL
#' 'kmer' -> one parameter necessary to set the length of the subsequence, defaults to 3
#' 'blosum' -> two parameters necessary: k ( The number of selected scales (i.e. the first k scales) derived by the substitution matrix. This can be selected according to the printed relative importance values.)
#' and lag (The lag parameter. Must be less than the amino acids.). They default to (5, 7).
#' 'topoPCA' -> three parameters necessary: index (Integer vector. Specify which molecular descriptors to select from the topological descriptors), pc (Integer. Number of principal components. Must be no greater than the number of amino acid properties provided.)
#' and lag(The lag parameter. Must be less than the amino acids.). They default to (c(1:78),5,7).
#' @param which.label String. The name of the column in VDJ which will be used as a label in a chosen model later. If missing, no label will be appended to the encoded features.
#' @param problem String ("classification" or "regression"). Whether the return matrix will be used in a classification problem or a regression one. Defaults to "classification".
#' @param verbose.classes Boolean. Whether to display information on the distribution of samples between classes. Defaults to TRUE. 
#' For this parameter to be set to TRUE, classification must all be set to TRUE (default).
#' @param platypus.version This function works with "v3" only, there is no need to set this parameter.
#' @return A dataframe containing the encoded features and its label, each row corresponding to a different cell.
#' The encodings are ordered as they have been entered in the 'which.features' parameter. 
#' The label can be found in the last column of the dataframe returned. 
#' @export
#' @examples
#' \dontrun{
#' To return the encoded 'VDJ_cdr3s_nt' sequences using 3mer encoding for nt sequences.
#' Attaching the "GP33_binder" label to be used in downstream ML models.
#' 
#' features_VDJ_GP33_binder <- PlatypusML_feature_extraction_VDJ(VGM = VGM, 
#' which.features = c("VDJ_cdr3s_nt"), 
#' which.encoding = c("kmer"), 
#' parameters.encoding.nt = c(3), 
#' which.label = "GP33_binder")
#'}


PlatypusML_feature_extraction_VDJ <- function(VGM,
                                          which.features,
                                          which.encoding,
                                          encoding.level,
                                          unique.sequence,
                                          parameters.encoding.nt,
                                          parameters.encoding.aa,
                                          which.label,
                                          problem,
                                          verbose.classes, 
                                          platypus.version){
  

  
  # checking for missing inputs and set defaults 
  if(missing(VGM)) stop("Please provide VDJ input for this function")
  if(missing(which.features)) stop("Please provide at  least one feature input for this function")
  if(missing(which.encoding)){
    which.encoding <- c("onehot", "onehot")
    print("which.encoding parameter was set to onehot")}
  if(length (which.encoding) != length (which.features)) stop("The input for which.ecoding has an incorrect size. A vector of the same size as which.features needs to be passed.")

  
  if (missing(encoding.level))  encoding.level <- "cell"
  if (encoding.level == "unique.sequence" & missing(unique.sequence))  unique.sequence <- "VDJ_cdr3s_aa"
  
  if (missing(parameters.encoding.nt)) parameters.encoding.nt<- NaN
  if (missing(parameters.encoding.aa)) parameters.encoding.aa<- NaN
  if (missing(which.label)) which.label <- NaN
  if (missing(problem)) problem <- "classification"
  if (missing(verbose.classes)) verbose.classes <- TRUE
  if (missing(platypus.version)) platypus.version <- "v3"
  
  ########################################################
  #extracting relevant matrix 
  VDJ <- VGM[[1]]
  
  ############### filtering ##############################
  
  #selecting the subset of data which has available labels 
  #in the case of classification, get rid of "nd" labels
  if (!is.na(which.label)){
    if (problem == "classification"  & nrow(VDJ[which(VDJ[which.label] =="nd"),][which.label])>=1)  {
    VDJ[which(VDJ[which.label] =="nd"),][which.label] <- NA
    }
  }
  
  #select only the part of the data that has available labels
  if (!is.na(which.label))  VDJ <- subset(VDJ, !is.na(VDJ[which.label]))
  
  #based on which encoding level you use filter the data 
  
  if (encoding.level == "unique.sequence"){
    VDJ <- VDJ %>% distinct (VDJ[unique.sequence], .keep_all = TRUE)
  }
  
  if (encoding.level == "clone"){
    VDJ$clonotype_distinct <- paste(VDJ$clonotype_id_10x, VDJ$sample_id)
    VDJ <- VDJ %>% distinct(VDJ$clonotype_distinct, .keep_all = TRUE)
  }
  
  #if the encoding level is "cell" no selection is required and the encoding will be done for all of the sequences 
  
  

  ############### define helper functions ################
  
  onehot <- function(sequences, sequence.type){
    aa <-  c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    nt <- c("A", "C", "G", "T")
    
    if (sequence.type == 'aa'){
      ref.seq <- aa
      ref.size <- 20
    }  else if (sequence.type == 'nt') {
      ref.seq <- nt
      ref.size <- 4
    }  
    #splitting the sequences 
    seq_split <- lapply(sequences, function(x)
      c(stringr::str_split_fixed(x, "", max(sapply(x, stringr::str_length)))))
    
    #identifying the maximum length of the sequence to be encoded 
    maximum = max(sapply(seq_split, function(x) length (x)))
    
    #pad the sequence to have the same length as the sequence with the maximum length 
    seq_pad <- lapply(seq_split, function(xx){
      if (maximum > length(xx)) {
        return( c(xx, rep(0, times=maximum-length(xx))) )
      } else {
        return(xx)
      }
    })
    
    #one-hot encoding
    x<-matrix(0,length(seq_pad),maximum*ref.size)
    #for each sequence define a matrix to store the encoding
    for(i in 1:length(seq_pad)){
      temp<-matrix(0,ref.size,length(seq_pad[[i]]))
      
      #encode each base in the sequence
      for(j in 1:length(seq_pad[[i]])){
        vec <- rep(0, times=ref.size)
        vec[which(ref.seq == seq_pad[[i]][j])]<-1
        temp[,j]<-vec
      }
      #flatten it and added to the final matrix
      x[i,]<-as.numeric(temp)
    }
    return (data.matrix(x))
    
  }
  
  ################################################
  
  #allocate array for encoded features
  encoding <- NULL
  
  ######## sequence encoding (nt or aa) ########### 
  
  for (i in 1:length(which.features)){
    
    feature <- which.features[i]
    which.encoding.feature <- which.encoding[i]
    
    #setting defaults
    if(is.null(which.encoding.feature)) which.encoding.feature <- "onehot"
    
    selected_feature <- VDJ[[feature]]
    
    #check if we are dealing with a nucleotide or aminoacid sequence 
    if (grepl('nt', feature, fixed = TRUE))  {
      seq.type <- 'nt'
      if(which.encoding.feature != "onehot" & which.encoding.feature != "kmer") stop("The wrong type of encoding was provided")
      if(missing(parameters.encoding.nt) &which.encoding.feature == "kmer") parameters.encoding.nt = c(3)
      parameters.encoding<-parameters.encoding.nt
    }
    else {
      seq.type <- 'aa' 
      if(which.encoding.feature != "onehot" & which.encoding.feature != "kmer" & which.encoding.feature != "blosum" &
         which.encoding.feature != "tc" & which.encoding.feature != "dc" & which.encoding.feature != "topoPCA") stop("The wrong type of encoding was provided")
      if (missing(parameters.encoding.aa)){
        
        if(which.encoding.feature == "kmer") parameters.encoding.aa = c(3)
        if (which.encoding.feature == "blosum") parameters.encoding.aa = c(5, 7)
        if (which.encoding.feature == "topoPCA") parameters.encoding.aa = list(c(1:78), 5, 7)
      }
      parameters.encoding <- parameters.encoding.aa
    }
    
    if (which.encoding.feature == "onehot"){
      encoded_feature <- onehot(selected_feature, sequence.type = seq.type)
    }
    
    if(stringr::str_detect(which.encoding.feature, "kmer")){
      if (parameters.encoding<10){
        num <- parameters.encoding
      }
      else { #23 or 234 kmer encoding
        num <- 2
      }
      for (feature in which.features){
        VDJ <- VDJ[nchar(VDJ[[feature]])>=num,]
      }
      
      if (seq.type == "nt"){
        num <- parameters.encoding.nt
        #simple kmer encoding strategy
        if (parameters.encoding<10){
          encoded_feature <- selected_feature %>%
            lapply(function(x) as.data.frame(kmer::kcount(ape::as.DNAbin(Biostrings::DNAStringSet(x)), k = num))) %>%
            bind_rows %>% data.matrix
          #normalize (to account for differences in sequence length)
          encoded_feature <- encoded_feature/rowSums(encoded_feature)
        }
        # concatenating 2mers 3mers encoding strategy
        if (parameters.encoding==23){
          twomers <- selected_feature %>%
            lapply(function(x) as.data.frame(kmer::kcount(ape::as.DNAbin(Biostrings::DNAStringSet(x)), k = 2))) %>%
            bind_rows %>% data.matrix
          #normalize (to account for differences in sequence length)
          twomers <- twomers/rowSums(twomers)
          threemers <- selected_feature %>%
            lapply(function(x) as.data.frame(kmer::kcount(ape::as.DNAbin(Biostrings::DNAStringSet(x)), k = 3))) %>%
            bind_rows %>% data.matrix
          #normalize (to account for differences in sequence length)
          threemers <- threemers/rowSums(threemers)
          encoded_feature <- cbind(twomers, threemers)
        }
        # concatenating 2mers 3mers 4mers encoding strategy
        if (parameters.encoding==234){
          twomers <- selected_feature %>%
            lapply(function(x) as.data.frame(kmer::kcount(ape::as.DNAbin(Biostrings::DNAStringSet(x)), k = 2))) %>%
            bind_rows %>% data.matrix
          #normalize (to account for differences in sequence length)
          twomers <- twomers/rowSums(twomers)
          threemers <- selected_feature %>%
            lapply(function(x) as.data.frame(kmer::kcount(ape::as.DNAbin(Biostrings::DNAStringSet(x)), k = 3))) %>%
            bind_rows %>% data.matrix
          #normalize (to account for differences in sequence length)
          threemers <- threemers/rowSums(threemers)
          fourmers <- selected_feature %>%
            lapply(function(x) as.data.frame(kmer::kcount(ape::as.DNAbin(Biostrings::DNAStringSet(x)), k = 4))) %>%
            bind_rows %>% data.matrix
          #normalize (to account for differences in sequence length)
          fourmers <- fourmers/rowSums(fourmers)
          encoded_feature <- cbind(twomers, threemers)
          encoded_feature <- cbind (encoded_feature, fourmers)
        }
      }
      
      if (seq.type == "aa"){
        num <- parameters.encoding.aa
        #simple kmer encoding strategy
        if (parameters.encoding<10){
          encoded_feature <- selected_feature %>%
            lapply(function(x) as.data.frame(kmer::kcount(ape::as.AAbin(Biostrings::AAStringSet(x)), k = num))) %>%
            bind_rows %>% data.matrix
          #normalize (to account for differences in sequence length)
          encoded_feature <- encoded_feature/rowSums(encoded_feature)
        }
        # concatenating 2mers 3mers encoding strategy
        if (parameters.encoding==23){
          twomers <- selected_feature %>%
            lapply(function(x) as.data.frame(kmer::kcount(ape::as.AAbin(Biostrings::AAStringSet(x)), k = 2))) %>%
            bind_rows %>% data.matrix
          #normalize (to account for differences in sequence length)
          twomers <- twomers/rowSums(twomers)
          threemers <- selected_feature %>%
            lapply(function(x) as.data.frame(kmer::kcount(ape::as.AAbin(Biostrings::AAStringSet(x)), k = 3))) %>%
            bind_rows %>% data.matrix
          #normalize (to account for differences in sequence length)
          threemers <- threemers/rowSums(threemers)
          fourmers <- selected_feature %>%
            lapply(function(x) as.data.frame(kmer::kcount(ape::as.AAbin(Biostrings::AAStringSet(x)), k = 4))) %>%
            bind_rows %>% data.matrix
          #normalize (to account for differences in sequence length)
          fourmers <- fourmers/rowSums(fourmers)
          encoded_feature <- cbind(twomers, threemers)
          encoded_feature <- cbind (encoded_feature, fourmers)
        }
      }
    }
    
    if(which.encoding.feature == "blosum"){
      k <- parameters.encoding[[1]]
      lag <- parameters.encoding[[2]]
      #used k and lag parameters as set by Damian
      encoded_feature <- selected_feature %>% lapply(function(x)as.data.frame(purrr::map(x, protr::extractBLOSUM, k=k, lag = lag )))%>% 
        bind_cols %>% t %>% matrix(nrow=length(selected_feature))
      
    }
    
    if(which.encoding.feature == "dc"){
      encoded_feature <- selected_feature %>% lapply(function(x)as.data.frame(purrr::map(x, protr::extractDC ))) %>% bind_cols %>% t %>% matrix(nrow=length(selected_feature))
      
    }
    
    if(which.encoding.feature == "tc"){
      encoded_feature <- selected_feature %>% lapply(function(x)as.data.frame(purrr::map(x, protr::extractTC ))) %>% bind_cols %>% t%>% matrix(nrow=length(selected_feature))
      
    }
    
    if(which.encoding.feature == "topoPCA"){
      index <- parameters.encoding[[1]]
      pc <- parameters.encoding[[2]]
      lag <- parameters.encoding[[3]]
      encoded_feature <- selected_feature %>% lapply(function(x)as.data.frame(purrr::map(x, protr::extractDescScales, propmat = "AATopo", index = index, pc = pc , lag = lag))) %>%
        bind_cols %>% t %>% matrix(nrow=length(selected_feature))
      
    }
    
    encoding <- cbind(encoding, encoded_feature)
  }
  
  #add label
  colnames(encoding) <- NULL #due to duplicate colnames when using kmer
  
  if (!is.na(which.label)){
    encoding <- cbind(encoding, VDJ[which.label])
  }

  #classes balance summary
  if (problem == "classification" & !is.na(which.label) & verbose.classes ){
    labels <- VDJ[which.label]
    classes <-  unique(labels)
    
    for (i in 1:nrow(classes)){
      count_class <- length(which(labels == classes[i,]))
      perc_class <- count_class/nrow(labels)*100
      print(sprintf("Class %s contains %s samples, representing %s%% of the total samples.", classes[i,], count_class, perc_class ))
    }
    
  }
  
  ################################################
  
  #return encoded extracted features 
  return (data.frame(encoding))
  
}