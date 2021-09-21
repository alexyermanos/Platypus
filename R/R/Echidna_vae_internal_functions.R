#onehot encoding-----
.onehot_encoding<-function(data,base){
  n<-stringr::str_length(data)
  pad<-c()
  l<-max(n)
  for(i in 1:length(n)) {
    if(n[i]<l){
      pad[i]<-paste0(rep(".",l-n[i]),collapse = "")
    } else pad[i]<-""
  }

  seq<-paste0(data,pad)
  data<-stringr::str_split(seq,pattern="")

  #one-hot encoding
  x<-matrix(0,length(data),l*4)
  for(i in 1:length(data)){
    temp<-matrix(0,4,length(data[[i]]))
    for(j in 1:length(data[[i]])){

      if(data[[i]][j]==base[1]) vec<-c(1,0,0,0)
      if(data[[i]][j]==base[2]) vec<-c(0,1,0,0)
      if(data[[i]][j]==base[3]) vec<-c(0,0,1,0)
      if(data[[i]][j]==base[4]) vec<-c(0,0,0,1)
      if(!(data[[i]][j] %in% base)) vec<-c(0,0,0,0)
      temp[,j]<-vec
    }
    x[i,]<-as.numeric(temp)
  }
  return(x)
}

#----Argmax----
.Arg_base<-function(matrix,null.threshold,base){
  seq<-c()
  for(i in 1:dim(matrix)[2]){
    if(max(matrix[,i])>=null.threshold){
      seq[i]<-base[which(matrix[,i]== max(matrix[,i]))][1]
    }
    if(mean(matrix[,i])<null.threshold){
      seq[i]<-""
    }
  }
  return(paste0(seq,collapse=""))
}

#----Arg one hot----
.Arg_onehot<-function(matrix,null.threshold){
  seq<-matrix

  for(i in 1:dim(matrix)[2]){
    vec<-c(0,0,0,0)
    if(max(matrix[,i])>=null.threshold){
      vec[which(matrix[,i]== max(matrix[,i]))]<-1
      seq[,i]<-vec
    }
    if(mean(matrix[,i])<null.threshold){
      seq[,i]<-vec
    }
  }
  return(seq)
}
