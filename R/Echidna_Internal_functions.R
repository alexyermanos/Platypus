#VDJ-----
.VDJ_RECOMBIN_FUNCTION <- function(v_seq, d_seq, j_seq,
                                   method, chain.type, species,
                                   vdj.insertion.mean,
                                   vdj.insertion.stdv){
  vdj_length_prob<-vdj_length_prob
  base_array <- c("a", "t", "g", "c")
  base_probs <- c(.25,.25,.25,.25)
  ## need to recombine V and D and D and J
  if (method=="naive" && chain.type=="heavy"){
    base_array <- c("a", "t", "g", "c")
    base_probs <- c(.25,.25,.25,.25)
    vdj_options <- c("none", "deletion", "insertion")
    vdj_options_prob <- c(1/3,1/3,1/3)
    recomb_decision <- sample(x = vdj_options, 2, replace=TRUE, prob = vdj_options_prob)

    ############ starts the combination between V and D
    if(recomb_decision[1]=="none"){ # no change
      v_d <- paste(v_seq, d_seq, sep="")
    }
    else if(recomb_decision[1]=="deletion"){
      # v_seq_new loses 3' end and d_seq_new1 loses up to 5' bases
      v_seq_new <- substring(v_seq,first=1, last=nchar(v_seq)-sample(x=vdj_length_prob[["v3_deletion"]]$length,1,replace=TRUE, prob=vdj_length_prob[["v3_deletion"]]$probability))
      d_seq_new5 <- substring(d_seq, first=sample(x=vdj_length_prob[["d5_deletion"]]$length+1,1,replace=TRUE,prob=vdj_length_prob[["d5_deletion"]]$probability), last=nchar(d_seq))
      v_d <- paste(v_seq_new, d_seq_new5, sep="")
    }
    else if(recomb_decision[1]=="insertion"){
      insertion_length <- sample(x=vdj_length_prob[["vd_insertion"]]$length,1,replace=TRUE, prob=vdj_length_prob[["vd_insertion"]]$probability)
      insertion_array <- sample(x=base_array, insertion_length, replace=TRUE, base_probs)
      insertion_string <- paste(insertion_array, collapse='')
      v_d <- paste(v_seq, insertion_string, d_seq, sep="")

    } ### END recombination between v and d regions and start D and J
    if(recomb_decision[2]=="none"){
      vdj <- paste(v_d, j_seq, sep="")
    }
    else if(recomb_decision[2]=="deletion"){
      ### need to cut off the end of v_d and cut off the 5' of J
      v_d_new <- substring(v_d,first=1, last=nchar(v_d)-sample(x=vdj_length_prob[["d3_deletion"]]$length,1,replace=TRUE, prob=vdj_length_prob[["d3_deletion"]]$probability))
      j_new <- substring(j_seq, first=sample(x=vdj_length_prob[["j5_deletion"]]$length,1,replace=TRUE,prob=vdj_length_prob[["j5_deletion"]]$probability), last=nchar(j_seq))
      vdj <- paste(v_d_new, j_new, sep="")
    }
    else if(recomb_decision[2]=="insertion"){
      dj_insertion_length <- sample(x=vdj_length_prob[["dj_insertion"]]$length,1,replace=TRUE, prob=vdj_length_prob[["dj_insertion"]]$probability)
      dj_insertion_array <- sample(x=base_array, dj_insertion_length, replace=TRUE, base_probs)
      dj_insertion_string <- paste(dj_insertion_array, collapse='')
      vdj <- paste(v_d, dj_insertion_string, j_seq, sep="")
    } ### ends the second recom_decision if-else


    return(vdj)
  }
  else if (chain.type=="light" && method=="naive"){
    base_array <- c("a", "t", "g", "c")
    base_probs <- c(.25,.25,.25,.25)
    vdj_options <- c("none", "deletion", "insertion")
    vdj_options_prob <- c(1/3,1/3,1/3)
    ## changed to only one event for light chain
    recomb_decision <- sample(x = vdj_options, 1, replace=TRUE, prob = vdj_options_prob)

    ############ starts the combination between V and D
    if(recomb_decision[1]=="none"){ # no change
      v_j <- paste(v_seq, j_seq, sep="")
    }
    else if(recomb_decision[1]=="deletion"){
      # v_seq_new loses 3' end and d_seq_new1 loses up to 5' bases
      v_seq_new <- substring(v_seq,first=1, last=nchar(v_seq)-sample(x=vdj_length_prob[["v3_deletion"]]$length,1,replace=TRUE, prob=vdj_length_prob[["v3_deletion"]]$probability))
      j_seq_new5 <- substring(j_seq, first=sample(x=vdj_length_prob[["j5_deletion"]]$length+1,1,replace=TRUE,prob=vdj_length_prob[["j5_deletion"]]$probability), last=nchar(j_seq))
      v_j <- paste(v_seq_new, j_seq_new5, sep="")
    }
    else if(recomb_decision[1]=="insertion"){
      insertion_length <- sample(x=vdj_length_prob[["vj_insertion"]]$length,1,replace=TRUE, prob=vdj_length_prob[["vj_insertion"]]$probability)
      insertion_array <- sample(x=base_array, insertion_length, replace=TRUE, base_probs)
      insertion_string <- paste(insertion_array, collapse='')
      v_j <- paste(v_seq, insertion_string, j_seq, sep="")

    } ### END recombination between v and j regions light chain

    return(v_j)
  }
  ## use data driven method for heavy chain with new insertion probabilities
  else if (method=="data" && chain.type=="heavy"){
    base_array <- c("a", "t", "g", "c")
    base_probs <- c(.25,.25,.25,.25)
    vdj_options <- c("none", "deletion", "insertion")
    vdj_options_prob <- c(1/3,1/3,1/3)
    recomb_decision <- sample(x = vdj_options, 2, replace=TRUE, prob = vdj_options_prob)

    ############ starts the combination between V and D
    if(recomb_decision[1]=="none"){ # no change
      v_d <- paste(v_seq, d_seq, sep="")
    }
    else if(recomb_decision[1]=="deletion"){
      # v_seq_new loses 3' end and d_seq_new1 loses up to 5' bases
      v_seq_new <- substring(v_seq,first=1, last=nchar(v_seq)-sample(x=c(0,1,2,3,4,5),1,replace=TRUE, prob=c(0.3,0.2,0.2,0.1,0.1,0.1)))
      d_seq_new5 <- substring(d_seq, first=sample(x=c(1,2,3,4,5,6),1,replace=TRUE,prob=c(0.3,0.2,0.2,0.1,0.1,0.1)), last=nchar(d_seq))
      v_d <- paste(v_seq_new, d_seq_new5, sep="")
    }
    else if(recomb_decision[1]=="insertion"){
      if(species=="mus" && vdj.insertion.mean=="default"){
        insertion_length <- stats::rnorm(n=1,mean=4.0,sd=vdj.insertion.stdv)
      }
      else if(species=="hum" && vdj.insertion.mean=="default"){ ### need to update this still
        insertion_length <- sample(x=c(1,2,3,4,5),1,replace=TRUE, prob=c(0.2,0.2,0.2,0.2,0.2))
      }
      else{
        insertion_length <- as.integer(abs(stats::rnorm(n=1,mean=vdj.insertion.mean,sd=vdj.insertion.stdv)))
      }
      #insertion_length <- sample(x=c(1,2,3,4,5),1,replace=TRUE, prob=c(0.2,0.2,0.2,0.2,0.2))
      insertion_array <- sample(x=base_array, insertion_length, replace=TRUE, base_probs)
      insertion_string <- paste(insertion_array, collapse='')
      v_d <- paste(v_seq, insertion_string, d_seq, sep="")

    } ### END recombination between v and d regions and start D and J
    if(recomb_decision[2]=="none"){
      vdj <- paste(v_d, j_seq, sep="")
    }
    else if(recomb_decision[2]=="deletion"){
      ### need to cut off the end of v_d and cut off the 5' of J
      v_d_new <- substring(v_d,first=1, last=nchar(v_d)-sample(x=c(1,2,3,4,5),1,replace=TRUE, prob=c(0.2,0.2,0.2,0.2,0.2)))
      j_new <- substring(j_seq, first=sample(x=c(1,2,3,4,5),1,replace=TRUE,prob=c(0.2,0.2,0.2,0.2,0.2)), last=nchar(j_seq))
      vdj <- paste(v_d_new, j_new, sep="")
    }
    else if(recomb_decision[2]=="insertion"){
      if(species=="mus" && vdj.insertion.mean=="default"){
        dj_insertion_length <- as.integer(abs(stats::rnorm(n=1,mean=2.9,sd=vdj.insertion.stdv)))
      }
      else if(species=="hum" && vdj.insertion.mean=="default"){ ### need to update this still for human dj
        dj_insertion_length <- sample(x=c(1,2,3,4,5),1,replace=TRUE, prob=c(0.2,0.2,0.2,0.2,0.2))
      }
      else{
        dj_insertion_length <- as.integer(abs(stats::rnorm(n=1,mean=vdj.insertion.mean,sd=vdj.insertion.stdv)))
      }
      #dj_insertion_length <- sample(x=c(2,4,6,8,10),1,replace=TRUE, prob=c(.2,.2,.2,.2,.2))
      dj_insertion_array <- sample(x=base_array, dj_insertion_length, replace=TRUE, base_probs)
      dj_insertion_string <- paste(dj_insertion_array, collapse='')
      vdj <- paste(v_d, dj_insertion_string, j_seq, sep="")
    } ### ends the second recom_decision if-else


    return(vdj)
  }
  else if (chain.type=="light" && method=="data"){
    base_array <- c("a", "t", "g", "c")
    base_probs <- c(.25,.25,.25,.25)
    vdj_options <- c("none", "deletion", "insertion")
    vdj_options_prob <- c(1/3,1/3,1/3)
    ## changed to only one event for light chain
    recomb_decision <- sample(x = vdj_options, 1, replace=TRUE, prob = vdj_options_prob)

    ############ starts the combination between V and D
    if(recomb_decision[1]=="none"){ # no change
      v_j <- paste(v_seq, j_seq, sep="")
    }
    else if(recomb_decision[1]=="deletion"){
      # v_seq_new loses 3' end and d_seq_new1 loses up to 5' bases
      v_seq_new <- substring(v_seq,first=1, last=nchar(v_seq)-sample(x=c(0,1,2,3,4,5),1,replace=TRUE, prob=c(0.3,0.2,0.2,0.1,0.1,0.1)))
      j_seq_new5 <- substring(j_seq, first=sample(x=c(1,2,3,4,5,6),1,replace=TRUE,prob=c(0.3,0.2,0.2,0.1,0.1,0.1)), last=nchar(j_seq))
      v_j <- paste(v_seq_new, j_seq_new5, sep="")
    }
    else if(recomb_decision[1]=="insertion"){
      if(vdj.insertion.mean=="default") insertion_length <- new_SHM_prob <- as.integer(abs(stats::rnorm(n=1,mean=4,sd=vdj.insertion.stdv)))
      else{
        insertion_length <- as.integer(abs(stats::rnorm(n=1,mean=vdj.insertion.mean,sd=vdj.insertion.stdv)))
      }
      #insertion_length <- sample(x=c(1,2,3,4,5),1,replace=TRUE, prob=c(0.2,0.2,0.2,0.2,0.2))
      insertion_array <- sample(x=base_array, insertion_length, replace=TRUE, base_probs)
      insertion_string <- paste(insertion_array, collapse='')
      v_j <- paste(v_seq, insertion_string, j_seq, sep="")

    } ### END recombination between v and j regions light chain

    return(v_j)
  }


}

#SHM(adapted)-----------
.SHM_FUNCTION_SEQUENCE4 <- function(vdj_seq, mut_param,v_seq,SHM.nuc.prob){

  ctrl <- 0
  if(mut_param=="naive" || mut_param=="all" ||mut_param=="poisson"){
    holding_mut <- sample(x=c(0,1), nchar(vdj_seq), replace=TRUE,
                          c(SHM.nuc.prob[1], 1-SHM.nuc.prob[1]))
    for (i in 1:nchar(vdj_seq)){
      if(holding_mut[i]==0){
        holding_char <- substr(vdj_seq, i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(rep(1/3,3)))
        #substr(vdj_seq, i, i) <- sample(x=c("A", "T", "G","C"), 1, replace=TRUE, c(rep(.25,4)))
        ctrl <- ctrl + 1
      }
    }
  }
  if(mut_param=="data" || mut_param=="all"){
    index_CDR_start <- nchar(v_seq)-15
    index_CDR_stop <- nchar(vdj_seq)
    CDR_length <- index_CDR_stop - index_CDR_start
    if(mut_param=="data") CDR_prob <- SHM.nuc.prob
    else if(mut_param=="all") CDR_prob <- SHM.nuc.prob[2]
    no_CDR_prob <- 1-CDR_prob
    #CDR_mut <- sample(x=c(0,1), nchar(CDR_length), replace=TRUE, c(no_CDR_prob,CDR_prob))
    for(i in 81:114){
      CDR_mut <- sample(x=c(0,1), 1, replace=TRUE, c(no_CDR_prob,CDR_prob))
      if(CDR_mut==1){
        holding_char <- substr(vdj_seq,i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(15,70,15))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(15,15,70))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(70,15,15))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(15,15,75))
        #substr(vdj_seq,i,i) <- sample(x=c("A","T","G","C"),1,replace=TRUE, c(rep(.25,4)))
        ctrl <- ctrl + 1
      }
    }
    for(i in 168:195){
      CDR_mut <- sample(x=c(0,1), 1, replace=TRUE, c(no_CDR_prob,CDR_prob))
      if(CDR_mut==1){
        holding_char <- substr(vdj_seq,i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(15,70,15))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(15,15,70))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(70,15,15))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(15,15,75))
        #substr(vdj_seq,i,i) <- sample(x=c("A","T","G","C"),1,replace=TRUE, c(rep(.25,4)))
        ctrl <- ctrl + 1
      }
    }
    for(i in index_CDR_start:index_CDR_stop){
      CDR_mut <- sample(x=c(0,1), 1, replace=TRUE, c(no_CDR_prob,CDR_prob))
      if(CDR_mut==1){
        holding_char <- substr(vdj_seq,i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(15,70,15))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(15,15,70))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(70,15,15))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(15,15,75))
        #substr(vdj_seq,i,i) <- sample(x=c("A","T","G","C"),1,replace=TRUE, c(rep(.25,4)))
        ctrl <- ctrl + 1
      }
    }
  }
  if(mut_param == "motif" || mut_param == "all"){
    hotspot_df<-hotspot_df
    random_spot <- hotspot_df[sample(nrow(hotspot_df)),]
    current_hot_spots <- 0
    if(mut_param=="motif") hot_spot_limit <- SHM.nuc.prob
    else if(mut_param=="all") hot_spot_limit <- SHM.nuc.prob[3]
    #hot_spot_limit <- SHM.nuc.prob
    for(i in 1:nrow(random_spot)){
      if(grepl(pattern = random_spot$pattern[i],x = vdj_seq)){
        current_hot_spots <- current_hot_spots+1
      }
      if(current_hot_spots>hot_spot_limit) break
      vdj_seq <- base::sub(random_spot$pattern[i],replacement = paste(substring(random_spot$pattern[i],first=1,last=2),sample(x = c("a","c","g","t"),1,replace=TRUE,prob = c(random_spot$toA[i], random_spot$toC[i], random_spot$toG[i], random_spot$toT[i])),substring(random_spot$pattern[i],first=4,last=5),sep=""),x = vdj_seq)
      ctrl<-ctrl+1
    }
  }
  if(mut_param == "wrc" || mut_param == "all"){
    one_spot_df<-one_spot_df
    random_spot <- one_spot_df[sample(nrow(one_spot_df)),]
    current_hot_spots <- 0
    if(mut_param=="wrc") hot_spot_limit <- SHM.nuc.prob
    else if(mut_param=="all") hot_spot_limit <- SHM.nuc.prob[4]
    for(i in 1:nrow(random_spot)){
      if(grepl(pattern = random_spot$pattern[i],x = vdj_seq)){
        current_hot_spots <- current_hot_spots+1
      }
      if(current_hot_spots>hot_spot_limit) break
      vdj_seq <- base::sub(random_spot$pattern[i],replacement = paste(substring(random_spot$pattern[i],first=1,last=2),sample(x = c("a","c","g","t"),1,replace=TRUE,prob = c(random_spot$toA[i], random_spot$toC[i], random_spot$toG[i], random_spot$toT[i])),substring(random_spot$pattern[i],first=4,last=5),sep=""),x = vdj_seq)
      ctrl<-ctrl+1
      }
  }
  return(list(vdj_seq,ctrl))

}


# .CLASS.SWITCH <- function(current_isotype,current_c_gene, class_switch_prob){
#
#   is_switch <- sample(x=c(0,1), replace=TRUE,size = 1, prob=c(class_switch_prob,1- class_switch_prob))
#   switched_c_gene<-current_c_gene
#   if (is_switch == 0 && current_isotype < "5"){
#     class.number<-1:(5-current_isotype)
#     current_isotype <- current_isotype+ sample(x=class.number, 1, replace=TRUE)
#     if (current_isotype==2){switched_c_gene<-"IGHD"}
#     else if(current_isotype==1){switched_c_gene<-"IGHM"}
#     else if (current_isotype==3){switched_c_gene<-"IGHG"}
#     else if (current_isotype==4){switched_c_gene<-"IGHE"}
#     else if (current_isotype==5){cswitched_c_gene<-"IGHA"}
#   }
#   return(list(current_isotype,switched_c_gene))
# }

.GENERATE.BARCODE<-function(size){
  bar<-0
  a<-1
  ATGC<-c("A","T","G","C")
  barcode_duo<-c()
  while(bar< size){
    code<-paste(paste(sample(ATGC,size=16, replace = T),collapse = ""),"-1",sep = "")
    if(code %in% barcode_duo){
      next}
    barcode_duo[2*a-1]<-code
    barcode_duo[2*a]<-code
    barcode_duo<-stats::na.omit(barcode_duo)
    bar<-length(barcode_duo)/2
    a<-a+1
  }
  return(barcode_duo)
}

.ADD.BARCODE<-function(barcode.duo){
  ATGC<-c("A","T","G","C")
  l<-length(barcode.duo)
  #new barcode
  bar<-0
  while(bar< l+1){
    code<-paste(paste(sample(ATGC,size=16, replace = T),collapse = ""),"-1",sep = "")
    if(code %in% barcode.duo){
      next}
    barcode.duo[l+1]<-code
    barcode.duo[l+2]<-code
    barcode.duo<-stats::na.omit(barcode.duo)
    bar<-length(barcode.duo)
  }
  return(barcode.duo)
}


#-----write consensus need to be followed by consensus<-consensus[,-1]
.WRITE.CONSENSUS<-function(clonotype.id,sequence.combined,barcode.unique,clonotypes.dataframe){
  raw_consensus_id<-c()
  consensus_seq<-c()
  ##single format generate consensus
  clonotype_seq<-data.frame(clonotype.id,sequence.combined,barcode.unique)
  colnames(clonotype_seq)<-c("clonotype_id","seq_combi", "barcode_uniq")
  clonotype_seq<-data.frame(merge(clonotypes.dataframe,clonotype_seq,by="clonotype_id",all=T))
  clonotype_seq$linked<-paste(clonotype_seq$clonotype.id,clonotype_seq$seq_combi,sep="_")
  #clonotype_seq: clonotype id, seq, clonotypefreq, porportion,linked clonotype and seq
  clonotype_freq<-reshape2::dcast(clonotype_seq,linked~.,value.var = "barcode_uniq",length)
  #ranking seq in each clonotype, seek for majority
  colnames(clonotype_freq)<-c("linked","seq_freq")
  #add linked freq to clonotype_seq
  clonotype_seq<-merge(clonotype_freq,clonotype_seq,by="linked",all=T)
  #rank the majority freq on the top
  clonotype_seq<-clonotype_seq[order(clonotype_seq$seq_freq,decreasing = T),]

  consent <- subset(clonotype_seq,!duplicated(clonotype_seq$clonotype_id),select = c("clonotype_id","seq_combi"))
  #delete the duplicated, only keep the majority seq
  #split heavy chain and light chain
  consent<-cbind(consent$clonotype_id, as.data.frame(stringr::str_split_fixed(consent$seq_combi, ";|:", n = 4)))
  clonotype_chain<-c()
  for(j in 1:nrow(consent)){
    raw_consensus_id[2*j-1]<-paste0(consent[j,1],"_consensus_1")
    raw_consensus_id[2*j]<-paste0(consent[j,1],"_consensus_2")
    consensus_seq[2*j-1]<-consent[j,3]
    consensus_seq[2*j]<-consent[j,5]
    clonotype_chain[2*j-1]<-paste0(consent[j,1],"_1")
    clonotype_chain[2*j]<-paste0(consent[j,1],"_2")
  }
  consent<-data.frame(clonotype_chain,raw_consensus_id,consensus_seq)
  return(consent)
}
#class switch(new)-----

# .CLASS.SWITCH<-function(current_isotype, class_switch_prob){
#   item<-c(as.numeric(which(class_switch_prob[current_isotype,]!=0)),current_isotype)
#   if(length(item)!=1){
#     iso_prob<-as.numeric(class_switch_prob[current_isotype,which(class_switch_prob[current_isotype,]!=0)])
#     who_switch <- sample(x=item, replace=T, size = 1, prob=c(iso_prob,1-sum(iso_prob)))#
#   }
#   if(length(item)==1){
#     who_switch <- current_isotype
#   }
#   is_switch<-current_isotype!=who_switch
#
#   return(list(who_switch,is_switch))
# }

##Class switch combined with trans.switch 2021.March9

.TRANS.SWITCH<-function(trans.state, trans.switch.prob){
  item<-c(as.numeric(which(trans.switch.prob[trans.state,]!=0)),trans.state)
  if(length(item)!=1){
    trans_prob<-as.numeric(trans.switch.prob[trans.state,which(trans.switch.prob[trans.state,]!=0)])
    who_switch <- sample(x=item, replace=T, size = 1, prob=c(trans_prob,1-sum(trans_prob)))
  }
  if(length(item)==1){
    who_switch <- trans.state
  }
  is_switch<-trans.state!=who_switch
  return(list(who_switch,is_switch))
}
#when trans.switch.SHM.dependant==T |trans.switch.isotype.dependant==T, transcriptome state for sure switch, but to which type is still decided by the probability matrix.
.TRANS.SWITCH.DEPENDANT<-function(trans.state,trans.switch.prob){
  item<-as.numeric(which(trans.switch.prob[trans.state,]!=0))
  if(length(item)==0){
    who_switch<-trans.state
  }
  if(length(item)==1){
    who_switch<-item
  }
  if(length(item)>1){
    trans_prob<-as.numeric(trans.switch.prob[trans.state,which(trans.switch.prob[trans.state,]!=0)])
    who_switch <- sample(x=item, replace=T, size = 1, prob=trans_prob)
  }
  is_switch<-T
  return(list(who_switch,is_switch))
}



.GENERATE.COMBI<-function(heavy.light.seq, chain.name){
  combi<-c()
  for(i in 1:length(heavy.light.seq)*0.5){
    combi[i]<-paste0(chain.name[2*i-1],":",heavy.light.seq[2*i-1],";",chain.name[2*i],":",heavy.light.seq[2*i])
  }
  return(combi)
}

.CELL.DIVISION.LINEAR<-function(clonotype_id,j,cell.division.prob){

clonofreq_temp<-as.data.frame(table(clonotype_id))
colnames(clonofreq_temp)<-c("clonotype_id", "Freq")
clonofreq_temp<-clonofreq_temp[order(clonofreq_temp$Freq),]

if(nrow(clonofreq_temp)>1){
  diff<-(cell.division.prob[2]-cell.division.prob[1])/(nrow(clonofreq_temp)-1)
  cell.division.prob.j<-(which(clonofreq_temp$clonotype_id==clonotype_id[j])-1)*diff+cell.division.prob[1]

}
if(nrow(clonofreq_temp)==1){
  cell.division.prob.j<-(cell.division.prob[2]+cell.division.prob[1])/2
}
is_cell_division <-  sample(x=c(0,1), replace=TRUE,size = 1, prob=c(cell.division.prob.j,1- cell.division.prob.j))
return(is_cell_division)
}

.CELL.DIVISION.LINEAR.INVERSE<-function(clonotype_id,j,cell.division.prob){

  clonofreq_temp<-as.data.frame(table(clonotype_id))
  colnames(clonofreq_temp)<-c("clonotype_id", "Freq")
  clonofreq_temp<-clonofreq_temp[order(clonofreq_temp$Freq,decreasing = T),]

  if(nrow(clonofreq_temp)>1){
    diff<-(cell.division.prob[2]-cell.division.prob[1])/(nrow(clonofreq_temp)-1)
    cell.division.prob.j<-(which(clonofreq_temp$clonotype_id==clonotype_id[j])-1)*diff+cell.division.prob[1]
  }
  if(nrow(clonofreq_temp)==1){
    cell.division.prob.j<-(cell.division.prob[2]+cell.division.prob[1])/2
  }
  is_cell_division <-  sample(x=c(0,1), replace=TRUE,size = 1, prob=c(cell.division.prob.j,1- cell.division.prob.j))
  return(is_cell_division)
}

.RM.EMPTY.NODE<-function(No,size0){
  if(length(size0)>0){
    size0<-rev(size0)
    empty_del<-as.data.frame(t(matrix(No,2)))
    for (e in size0){
      mom<-empty_del[,1][which(empty_del[,2]==e)]
      empty_del<-empty_del[-which(empty_del[,2]==e),]
      empty_del[,1][which(empty_del[,1]==e)]<-mom
      empty_del[empty_del>e]<-empty_del[empty_del>e]-1
    }
    No<-as.vector(t(empty_del))
  }

  return(No)
}




