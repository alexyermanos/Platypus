#' Plot the top abundant clonal frequencies in a barplot with ggplot2
#' @title Plot clonal frequency barplot of the outout simulated data
#' @param clonotypes The clonotypes dataframe, which is the second element in the simulation output list.
#' @param top.n The top n abundant clones to be shown in the plot. If missing, all clones will be shown.
#' @param y.limit The upper limit for y axis in the plot.
#' @return top abundant clonal frequencies in a barplot with ggplot2
#' @export
clonofreq<-function(clonotypes, top.n, y.limit){
  Clonotypes<-frequency <- NULL
  #require(ggplot2,dplyr)
  clonotypes<-clonotypes[order(clonotypes[,2],decreasing = T),]
  if(!(missing(top.n))) {
    clonotypes<-clonotypes[1:top.n,]
  }
  clonotypes$Clonotypes<-factor(c(1:nrow(clonotypes)))
  p<- ggplot2::ggplot( clonotypes, ggplot2::aes(x=Clonotypes, y=frequency)) +
    ggplot2::xlab("Clonotypes")+
    ggplot2::ylab("Cells")+
    ggplot2::geom_bar(stat="identity")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_x_discrete(expand = c(0,0))

    if(!missing(y.limit)){

    p<- ggplot2::ggplot( clonotypes, ggplot2::aes(x=Clonotypes, y=frequency)) +
    ggplot2::xlab("Clonotypes")+
    ggplot2::ylab("Cells")+
    ggplot2::geom_bar(stat="identity")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_x_discrete(expand = c(0,0))+
    ggplot2::scale_y_continuous(expand = c(0,0),limit = c(0, y.limit))

    }


  return(p)
}

#' Plot a stacked barplot for clonotype counts grouped by isotype.
#' @title Get information about the clonotype counts grouped by isotype.
#' @param all.contig.annotations The output dataframe all_contig_annotation from function simulate.repertoire.
#' @param top.n The top n abundant clones to be shown in the plot. If missing, all clones will be shown.
#' @param colors A named character vector of colors, the names are the isotypes. If missing, the default has 11 colors coresponding to the default isotype names.
#' @param y.limit The upper limit for y axis in the plot.
#' @return a stacked barplot for clonotype counts grouped by isotype
#' @export
clonofreq.isotype.plot<-function(all.contig.annotations,top.n,y.limit,colors){
  #require(dplyr,ggplot2)
  Clonotypes<-c_gene<-count<-NULL
  if(missing(colors))  {
    colors<-c("#377EB8","#4DAF4A",	"#984EA3",	"#FF7F00",	"#FFFF33",	"#A65628",	"#F781BF",	"#999999","#E41A1C","#984EA3",	"#FFFF33")
    names(colors)<-c("IGHM","IGHD","IGHG1","IGHG2A","IGHG2B","IGHG2C","IGHG3","IGHE",	"IGHA",	"IGHG2",	"IGHG4")
  }
  heavy_df<-all.contig.annotations[grep(all.contig.annotations$c_gene,pattern = "IGH"),]
  clonotypes<-as.data.frame(table(heavy_df$raw_clonotype_id))
  colnames(clonotypes)[1]<-"raw_clonotype_id"
  clonotypes<-clonotypes[order(clonotypes$Freq,decreasing = T),]
  clonotypes$Clonotypes<-factor(c(1:nrow(clonotypes)))
  heavy_df<-merge(heavy_df,clonotypes,all = T)
    if(!missing(top.n)){
    clonotypes<-clonotypes[1:top.n,]
    heavy_df<-heavy_df[heavy_df$raw_clonotype_id %in% clonotypes$raw_clonotype_id,]
    }


  iso<-dplyr::summarise(dplyr::group_by(heavy_df,Clonotypes,c_gene),count = dplyr::n())
  p<-ggplot2::ggplot(iso,ggplot2::aes(x=Clonotypes, y=count, fill = c_gene))+
    ggplot2::ylab("Cells")+
    ggplot2::xlab("Clonotypes")+
    ggplot2::geom_bar(stat="identity")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_x_discrete(expand = c(0,0))+
    ggplot2::scale_y_continuous(expand = c(0,0),limit = c(0, y.limit))+
    ggplot2::scale_fill_manual(values = colors,name = "Isotype")
  return(p)
}

#'Return 
#' @title Get information about the clonotype counts grouped by isotype.
#' @param all.contig.annotations The output dataframe all_contig_annotation from function simulate.repertoire.
#' @param top.n The top n abundant clones to be shown in the plot. If missing, all clones will be shown.
#' @return dataframes containing the top n abundant clonotypes and their frequency and isotype information for further processing.
#' @export
clonofreq.isotype.data<-function(all.contig.annotations,top.n){
  #require(dplyr)
  raw_clonotype_id<-c_gene<-NULL
  heavy_df<-all.contig.annotations[grep(all.contig.annotations$c_gene,pattern = "IGH"),]

  clonotypes<-as.data.frame(table(heavy_df$raw_clonotype_id))
  colnames(clonotypes)[1]<-"raw_clonotype_id"
  heavy_df<-merge(heavy_df,clonotypes,all = T)
  clonotypes<-clonotypes[order(clonotypes$Freq,decreasing = T),]

  if(!missing(top.n)){
    clonotypes<-clonotypes[1:top.n,]
    heavy_df<-heavy_df[heavy_df$raw_clonotype_id %in% clonotypes$raw_clonotype_id,]
  }

  iso<-dplyr::summarise(dplyr::group_by(heavy_df,raw_clonotype_id,c_gene),count = dplyr::n())
  return(iso)
}


#' Return a barplot of mean and standard error bar of certain value of each clone.
#' @param data A dataframe. Columns are different simulations, rows are the top clones. The first row is the top abundant clone.
#' @param y.lab A string specifies the y lable name of the barplot.
#' @param y.limit The upper limit for y axis in the plot.
#' @return a barplot of mean and standard error bar of certain value of each clone.
#' @export
get.barplot.errorbar<-function(data,y.lab,y.limit){
  sd<-NULL
  #require(ggplot2)
  Clones<-as.factor(c(1:nrow(data)))
  Mean<-rowMeans(data)
  se<-apply(data, 1,sd)/ncol(data)^0.5
  df<-data.frame(Clones,Mean,se)
  p<-ggplot2::ggplot(df,ggplot2::aes(x=Clones, y=Mean)) +
    ggplot2::ylab(y.lab)+
    ggplot2::xlab("Clonotypes")+
    ggplot2::geom_bar(stat="identity",  position=ggplot2::position_dodge())+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=Mean-se, ymax=Mean+se), width=.2,position=ggplot2::position_dodge(.9))+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_x_discrete(expand = c(0,0))+
    ggplot2::scale_y_continuous(expand = c(0,0),limit = c(0, y.limit))
  return(p)
}

#' Get the index of top ranking clones.
#' @param clonotypes The output "clonotypes" dataframe from simulation output.
#' @param top.n The top n abundant clones to be selected.
#' @return a vector of indexes of top ranking clones
#' @export
select.top<-function(clonotypes,top.n){
  #require(dplyr,stringr)
  clonotypes<-clonotypes[order(clonotypes$frequency,decreasing = T),][1:top.n,]
  selected<-as.numeric(stringr::str_replace(string=clonotypes$clonotype_id,pattern = "clonotype",replacement = ""))
  return(selected)
}

#' Get the number of unique variants in each clone in a vector and the barplot. The first item in the output is the vector representing the numbers of unique variants, the second item is the barplot.
#' @param igraph.index.attr The output "igraph.index.attr" list from simulation output.
#' @param clonotype.select The index of the clones to be shown. If missing, all clones will be included.
#' @param y.limit The upper limit for y axis in the plot.
#' @return the number of unique variants in each clone in a vector and the barplot. The first item in the output is the vector representing the numbers of unique variants, the second item is the barplot.
#' @export
get.n.node.plot<-function(igraph.index.attr,clonotype.select,y.limit){
 #require(ggplot2)
  Clonotypes<-unique_variants<-NULL
  n_node<-c()
  for(i in 1:length(igraph.index.attr)){
    if(length(igraph.index.attr[[i]])>1){
    n_node[i]<-nrow(unique(igraph.index.attr[[i]]))
    }else
      n_node[i]<-1
  }
  if(!(missing(clonotype.select))){
    n_node<-n_node[clonotype.select]
  }
  df<-data.frame(factor(c(1:length(n_node))),n_node)
  colnames(df)<-c("Clonotypes","unique_variants")
  p<- ggplot2::ggplot(df, ggplot2::aes(x=Clonotypes, y=unique_variants)) +
    ggplot2::xlab("Clonotypes")+
    ggplot2::ylab("Unique Variants")+
    ggplot2::geom_bar(stat="identity")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_x_discrete(expand = c(0,0))+
    ggplot2::scale_y_continuous(expand = c(0,0),limit = c(0, y.limit))
  return(p)
}

#' Get the number of unique variants in each clone in a vector. The output is the vector representing the numbers of unique variants.
#' @param data The output "igraph.index.attr" list from simulation output.
#' @param clonotype.select The index of the clones to be shown. If missing, all clones will be included.
#' @return the number of unique variants in each clone in a vector. The output is the vector representing the numbers of unique variants.
#' @export
get.n.node.data<-function(data,clonotype.select){
  n_node<-c()
  for(i in 1:length(data)){
    if(length(data[[i]])>1){
      n_node[i]<-nrow(unique(data[[i]]))
    }else
      n_node[i]<-1
  }
  if(!(missing(clonotype.select))){
    n_node<-n_node[clonotype.select]
  }
  df<-data.frame(paste0("clonotype",clonotype.select),n_node)
  colnames(df)<-c("clonotype_id","unique_variants")
  return(n_node)
}

#' Get the seurat object from simulated transciptome output.
#' @param data The output "transcriptome" dataframe from simulation output.
#' @return the seurat object from simulated transciptome output.
#' @export

get.elbow<-function(data){
  gex <- Seurat::CreateSeuratObject(counts = data)
  gex <- Seurat::NormalizeData(gex, normalization.method = "LogNormalize", scale.factor = 10000)
  gex <- Seurat::FindVariableFeatures(gex, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(gex)
  gex <- Seurat::ScaleData(gex,features = all.genes)
  gex <- Seurat::RunPCA(gex, features = Seurat::VariableFeatures(object = gex))
  return(gex)
}

#' Further process the seurat object from simulated transciptome output and make UMAP ready for plotting.
#' @param gex output from get.elbow function.
#' @param d dims argurment of in Seurat::FindNeighbors() and Seurat::RunUMAP
#' @param reso resolution argument in Seurat::FindClusters()
#' @return Further processed seurat object from simulated transciptome output with UMAP ready for plotting.
#' @export
get.umap<-function(gex,d,reso){
  #require(Seurat)
  gex <- Seurat::FindNeighbors(gex, dims = 1:d)
  gex <- Seurat::FindClusters(gex, resolution = reso)
  gex <- Seurat::RunUMAP(gex, dims = 1:d)
  return(gex)
}





#' Computing sequence distance according to the number of unmatched bases.
#' @param germline A string representing the germline sequence.
#' @param sequence  A string of the sequence to be compared, which has the same length as germline.
#' @return the number of unmatched bases in 2 sequences.
#' @export
get.seq.distance<-function(germline, sequence){

  a<-strsplit(germline,split = "")[[1]]
  b<-strsplit(sequence,split = "")[[1]]
  dis<-sum(a!=b)
  return(dis)
}

#' Get information about somatic hypermutation in the simulation. This function return a barplot showing the average mutation.
#' @param igraph.index.attr A list "igraph.index.attr" from the simulation output.
#' @param history  A dataframe "history" from the simulation output.
#' @param clonotype.select The selected clonotype index, can be the output of the function "select.top".
#' @param level Can be "node" or "cell". If "node", the function will return average mutation on unique variant level. Otherwise it will return on cell level.
#' @param y.limit The upper limit for y axis in the plot.
#' @return a barplot showing the average mutation per node (same heavy and light chain set) or per cell.
#' @export
get.avr.mut.plot<-function(igraph.index.attr,history,clonotype.select,level, y.limit){
  #require(stats,ggplot2)
  Clonotypes<-Mean<-Se<-NULL
  history<-as.data.frame(history)
  mut_by_node_mean<-c()
  mut_by_node_se<-c()
  mut_by_cell_mean<-c()
  mut_by_cell_se<-c()

  for (u in 1:length(clonotype.select)){

    dis_vec<-c()
    i<-clonotype.select[u]
    if ( length(igraph.index.attr[[i]])>1){
      igraph.index.attr[[i]]<-igraph.index.attr[[i]][-1,]
      rownames(igraph.index.attr[[i]])<-c(1:nrow(igraph.index.attr[[i]]))
      germ_seq_num<-igraph.index.attr[[i]][1,1]
      germ_seq<-history[which(history[,1]==germ_seq_num)[1],3]
      mut_seq_num<-igraph.index.attr[[i]][-1,1]

      for (j in 1:length(mut_seq_num)){
        mut_seq<-history[which(history[,1]==mut_seq_num[j])[1],3]
        dis_vec[j]<-get.seq.distance(germ_seq,mut_seq)
      }
      mut_by_node_mean[u]<-mean(dis_vec)
      mut_by_node_se[u]<-stats::sd(dis_vec)/length(dis_vec)^0.5
      mut_by_cell_mean[u]<-mean(rep(dis_vec,igraph.index.attr[[i]][["size"]][-1]))
      mut_by_cell_se[u]<-stats::sd(rep(dis_vec,igraph.index.attr[[i]][["size"]][-1]))/sum(igraph.index.attr[[i]][["size"]][-1])^0.5
    }

    if ( length(igraph.index.attr[[i]])==1){
      mut_by_node_mean[u]<-0
      mut_by_node_se[u]<-0
      mut_by_cell_mean[u]<-0
      mut_by_cell_se[u]<-0
    }
  }
  if(level=="node"){
  mut<-data.frame(clonotype.select,mut_by_node_mean,mut_by_node_se)
  }
  if(level=="cell"){
  mut<-data.frame(clonotype.select,mut_by_cell_mean,mut_by_cell_se)
  }
  colnames(mut)<-c("clonotype_id","Mean","Se")
  mut$Clonotypes<-factor(c(1:nrow(mut)))
  p<- ggplot2::ggplot(mut, ggplot2::aes(x=Clonotypes, y=Mean)) +
    ggplot2::xlab("Clonotypes")+
    ggplot2::ylab("Average mutation")+
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge())+
    ggplot2::geom_errorbar( ggplot2::aes(ymin=Mean-Se, ymax=Mean+Se), width=.2,
                            position= ggplot2::position_dodge(.9))+ ggplot2::theme_minimal()+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x =  ggplot2::element_text(angle = 45, hjust = 1))+
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_x_discrete(expand = c(0,0))+
    ggplot2::scale_y_continuous(expand = c(0,0),limit = c(0, y.limit))
  return(p)
}

#' Get information about somatic hypermutation in the simulation. This function return a barplot showing the average mutation.
#' @param igraph.index.attr A list "igraph.index.attr" from the simulation output.
#' @param history  A dataframe "history" from the simulation output.
#' @param clonotype.select The selected clonotype index, can be the output of the function "select.top".
#' @param level Can be "clone" or "cell". If "clone", the function will return average mutation on unique variant level. Otherwise it will return on cell level.
#' @return a bar plot showing the average mutation on clone or cell level.
#' @export
get.avr.mut.data<-function(igraph.index.attr,history,clonotype.select,level){
  #require(stats)
  history<-as.data.frame(history)
  mut_by_node_mean<-c()
  mut_by_node_se<-c()
  mut_by_cell_mean<-c()
  mut_by_cell_se<-c()

  for (u in 1:length(clonotype.select)){
    dis_vec<-c()
    i<-clonotype.select[u]

    if ( length(igraph.index.attr[[i]])>1){
      igraph.index.attr[[i]]<-igraph.index.attr[[i]][-1,]
      rownames(igraph.index.attr[[i]])<-c(1:nrow(igraph.index.attr[[i]]))

      germ_seq_num<-igraph.index.attr[[i]][1,1]
      germ_seq<-history[which(history[,1]==germ_seq_num)[1],3]
      mut_seq_num<-igraph.index.attr[[i]][-1,1]

      for (j in 1:length(mut_seq_num)){
        mut_seq<-history[which(history[,1]==mut_seq_num[j])[1],3]
        dis_vec[j]<-get.seq.distance(germ_seq,mut_seq)
      }
      mut_by_node_mean[u]<-mean(dis_vec)
      mut_by_node_se[u]<-stats::sd(dis_vec)/length(dis_vec)^0.5
      mut_by_cell_mean[u]<-mean(rep(dis_vec,igraph.index.attr[[i]][["size"]][-1]))
      mut_by_cell_se[u]<-stats::sd(rep(dis_vec,igraph.index.attr[[i]][["size"]][-1]))/sum(igraph.index.attr[[i]][["size"]][-1])^0.5
    }

    if ( length(igraph.index.attr[[i]])==1){
      mut_by_node_mean[u]<-0
      mut_by_node_se[u]<-0
      mut_by_cell_mean[u]<-0
      mut_by_cell_se[u]<-0
    }
  }
  if(level=="node"){
    mut<-data.frame(clonotype.select,mut_by_node_mean,mut_by_node_se)
  }
  if(level=="cell"){
    mut<-data.frame(clonotype.select,mut_by_cell_mean,mut_by_cell_se)
  }
  return(mut)
}

#' Get paired v gene heavy chain and light chain matrix on clonotype level. A v gene usage pheatmap can be obtain by p<-pheatmap::pheatmap(vgu_matrix,show_colnames= T, main = "V Gene Usage"), where the vgu_matrix is the output of this function.
#' @param all.contig.annotations The dataframe "all_contig_annotation" from simulation output.
#' @param level Can be "clone" or "cell". If "clone", the function will return paired v gene usage matrix on clonotype level. Otherwise it will return on cell level.
#' @return a paired v gene heavy chain and light chain matrix on clonotype level.
#' @export
get.vgu.matrix<-function(all.contig.annotations, level){
  #require(stringr,reshape2)
  h_df<-subset(all.contig.annotations, stringr::str_detect(pattern="contig_1",string=all.contig.annotations$raw_contig_id),select= c("barcode","raw_clonotype_id","v_gene"))
  l_df<-subset(all.contig.annotations, stringr::str_detect(pattern="contig_2",string=all.contig.annotations$raw_contig_id),select= c("barcode","v_gene"))
  colnames(h_df)<-stringr::str_replace(string = colnames(h_df),pattern = "v_gene",replacement = "v_h")
  colnames(l_df)<-stringr::str_replace(string = colnames(l_df),pattern = "v_gene",replacement = "v_l")
  df<-merge(h_df,l_df,by="barcode")
  if(level=="clone"){
    vgu_matrix<-reshape2::dcast(unique(df[,-1]),v_h~v_l,value.var = "raw_clonotype_id",length)
  }
  if(level=="cell"){
    vgu_matrix<-reshape2::dcast(df,v_h~v_l,value.var = "barcode",length)
  }
  rownames(vgu_matrix)<-vgu_matrix[,1]
  vgu_matrix<-as.matrix(vgu_matrix[,-1])
  return(vgu_matrix)
}

#' Plot a stacked barplot for clonotype counts grouped by transcriptome state(cell type).
#' @title Get information about the clonotype counts grouped by transcriptome state(cell type).
#' @param all.contig.annotations The output dataframe all_contig_annotation from function simulate.repertoire.
#' @param history The dataframe history from simulate output.
#' @param trans.names The names of cell types which are used in transcriptome.switch.prob argument in the simulation.
#' @param top.n The top n abundant clones to be shown in the plot. If missing, all clones will be shown.
#' @param y.limit The upper limit for y axis in the plot.
#' @param colors A named character vector of colors, the names are the isotypes. If missing, the default has 11 colors coresponding to the default isotype names.
#' @return a stacked barplot for clonotype counts grouped by transcriptome state(cell type).
#' @export
clonofreq.trans.plot<-function(all.contig.annotations,history,trans.names,top.n,y.limit,colors){
  #require(stringr,dplyr,ggplot2)
  Clonotypes<-trans_names<-count<-NULL
  if(missing(colors)) {
    colors<-c("#377EB8","#4DAF4A",	"#984EA3",	"#FF7F00","#FFFF33", "#A65628", "#F781BF", "#999999")[1:length(trans.names)]
    names(colors)<-trans.names
  }
  trans.names<-data.frame(c(1:length(trans.names)),trans.names)
  colnames(trans.names)<-c("trans_state_his","trans_names")
  heavy_df<-all.contig.annotations[grep(all.contig.annotations$c_gene,pattern = "IGH"),]
  history<-as.data.frame(history)
  colnames(history)[stringr::str_detect(pattern="barcode",string=colnames(history))]<-"barcode"

  clonotypes<-as.data.frame(table(heavy_df$raw_clonotype_id))
  colnames(clonotypes)[1]<-"raw_clonotype_id"
  clonotypes<-clonotypes[order(clonotypes$Freq,decreasing = T),]
  clonotypes$Clonotypes<-factor(c(1:nrow(clonotypes)))

  heavy_df<-merge(heavy_df,clonotypes,all = T)
  heavy_df<-merge(heavy_df, history, by="barcode",all.x = T)
  heavy_df<-merge(heavy_df, trans.names,all.x = T)
  if(!missing(top.n)){
    clonotypes<-clonotypes[1:top.n,]
    heavy_df<-heavy_df[heavy_df$raw_clonotype_id %in% clonotypes$raw_clonotype_id,]
  }
  iso<-dplyr::summarise(dplyr::group_by(heavy_df,Clonotypes,trans_names),count = dplyr::n())
  p<-ggplot2::ggplot(iso,ggplot2::aes(x=Clonotypes, y=count, fill = "trans_names"))+
    ggplot2::ylab("Cells")+
    ggplot2::xlab("Clonotypes")+
    ggplot2::geom_bar(stat="identity")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_x_discrete(expand = c(0,0))+
    ggplot2::scale_y_continuous(expand = c(0,0),limit = c(0, y.limit))+
    ggplot2::scale_fill_manual(values = colors,name="Cell Types")
  return(p)
}

#' Dataframe with clonotype counts grouped by transcriptome state(cell type).
#' @title Get information about the clonotype counts grouped by transcriptome state(cell type).
#' @param all.contig.annotations The output dataframe all_contig_annotation from function simulate.repertoire.
#' @param history The dataframe history from simulate output.
#' @param trans.names The names of cell types which are used in transcriptome.switch.prob argument in the simulation.
#' @param top.n The top n abundant clones to be shown in the plot. If missing, all clones will be shown.
#' @return a dataframe with clonotype counts grouped by transcriptome state(cell type).
#' @export
clonofreq.trans.data<-function(all.contig.annotations,history,trans.names,top.n){
  #require(stringr,dplyr)
  Clonotypes<-trans_names<-NULL
  trans.names<-data.frame(c(1:length(trans.names)),trans.names)
  colnames(trans.names)<-c("trans_state_his","trans_names")
  heavy_df<-all.contig.annotations[grep(all.contig.annotations$c_gene,pattern = "IGH"),]
  history<-as.data.frame(history)
  colnames(history)[stringr::str_detect(pattern="barcode",string=colnames(history))]<-"barcode"

  clonotypes<-as.data.frame(table(heavy_df$raw_clonotype_id))
  colnames(clonotypes)[1]<-"raw_clonotype_id"
  clonotypes<-clonotypes[order(clonotypes$Freq,decreasing = T),]
  clonotypes$Clonotypes<-factor(c(1:nrow(clonotypes)))

  heavy_df<-merge(heavy_df,clonotypes,all = T)
  heavy_df<-merge(heavy_df, history, by="barcode",all.x = T)
  heavy_df<-merge(heavy_df, trans.names,all.x = T)
  if(!missing(top.n)){
    clonotypes<-clonotypes[1:top.n,]
    heavy_df<-heavy_df[heavy_df$raw_clonotype_id %in% clonotypes$raw_clonotype_id,]
  }
  iso<-dplyr::summarise(dplyr::group_by(heavy_df,Clonotypes,trans_names),count = dplyr::n())
  return(iso)
}



#' Set idents for top abundant clones in Seurat object, get ready for highlight the top abundant clones in UMAP.
#' @param gex output from get.umap function.
#' @param all.contig.annotations The output dataframe all_contig_annotations from simulation.
#' @param top.n The top n abundant clones to be shown in the plot. If missing, all clones will be shown.
#' @return a Seurat object ready for highlight the top abundant clones in UMAP
#' @export

umap.top.highlight<-function(gex,all.contig.annotations,top.n){
  #require(dplyr,Seurat)
  clonotypes<-as.data.frame(table(all.contig.annotations$raw_clonotype_id))
  colnames(clonotypes)[1]<-"raw_clonotype_id"
  clonotypes<-clonotypes[order(clonotypes$Freq,decreasing = T),]
  if(!missing(top.n)){
    clonotypes<-clonotypes[1:top.n,]
  }
  top<-list()
  Seurat::Idents(gex)<-"Unselected"
  for (i in 1: nrow(clonotypes)){
    top[[i]]<-all.contig.annotations$barcode[which(all.contig.annotations$raw_clonotype_id==clonotypes$raw_clonotype_id[i])]
    gex<- Seurat::SetIdent(object = gex, cells = top[[i]], value =paste0("CloneRank",i))
  }
  return(gex)
}

#' Get clone network igraphs colored by seurat cluster id.
#' @param meta.data the meta.data dataframe from the Seurat object of the simulation. The object should be pre-processed and has cluster ids in the meta.data.
#' @param history The dataframe 'history' from the simulation output.
#' @param igraph.index The list 'igraph.index' from the simulation output.
#' @param empty.node If TRUE, there will be empty node in igraph. if FALSE, the empty node will be deleted.
#' @return a list of clone network igraphs colored by seurat cluster id
#' @export
cluster.id.igraph<-function(meta.data,history,igraph.index,empty.node){
  #require(reshape2,stats,igraph)
  colors<-colors
  igraph.index.attr<-list()
  igraph.index.jr<-list()
  igraph_list<-list()
  pie.values.list<-list()
  igraph_list_trans<-list()
  pie_values_trans_list<-list()
  #size.of.vertex
  size.of.vertex<-as.data.frame(stats::aggregate(barcode.history~seq.number,history,length,na.action = stats::na.pass))
  size.of.vertex1<-as.data.frame(stats::aggregate(barcode.history~seq.number,history,length))
  size.of.vertex<-merge(size.of.vertex1,size.of.vertex,all=T, by="seq.number")[,-3]
  size.of.vertex[!(size.of.vertex$seq.number %in% size.of.vertex1$seq.number),2]<-0
  colnames(size.of.vertex)<-c("seq.number","size")
  rm(size.of.vertex1)
  #meta.data get cluster id
  meta.data$barcode.history<-rownames(meta.data)
  meta.data<-data.frame(meta.data$barcode.history,meta.data$seurat_clusters)
  colnames(meta.data)<-c("barcode.history","seurat.clusters")
  history<-merge(history,meta.data,by="barcode.history",all=T)

  cluster.distribution<-reshape2::dcast(history, seq.number~seurat.clusters, value.var ="seurat.clusters",length)[,-1]
  cluster.distribution<-cluster.distribution[,1:ncol(cluster.distribution)-1]
  cluster.name<-colnames(cluster.distribution)
  cluster.distribution <-as.list( as.data.frame(t(cluster.distribution)))
  #write igraph list and assign attr
  for (i in 1:length(igraph.index)){

    if (length(igraph.index[[i]])>1){
      igraph.index.attr[[i]]<-data.frame(igraph.index[[i]])
      colnames(igraph.index.attr[[i]])<-"seq.number"

      igraph.index.jr[[i]]<-data.frame(igraph.index[[i]])
      colnames(igraph.index.jr[[i]])<-"seq.number"

      #delete duplicated for vertex attribute
      igraph.index.attr[[i]]<-subset(igraph.index.attr[[i]],!duplicated(igraph.index.attr[[i]]$seq.number))
      #size of vertex
      igraph.index.attr[[i]]<-merge(size.of.vertex,igraph.index.attr[[i]],by="seq.number",all.y=T)
      igraph.index.attr[[i]]$size[1]<-0
      igraph.index.attr[[i]]$adj_size<-(igraph.index.attr[[i]]$size*60/(range(igraph.index.attr[[i]]$size)[2]-range(igraph.index.attr[[i]]$size)[1]))^0.5+15 #visual size of vertex a number between 20 and 60
      #size of empty
      igraph.index.attr[[i]]$adj_size<-replace(igraph.index.attr[[i]]$adj_size, igraph.index.attr[[i]]$size==0, 10)

      #simplified index
      igraph.index.jr[[i]]$No<-match(x=igraph.index.jr[[i]]$seq.number,table=igraph.index.attr[[i]]$seq.number)

      #find empty
      size0<-which(igraph.index.attr[[i]]$size==0)

      No<- igraph.index.jr[[i]]$No

      #option: no empty node
      if(empty.node==F&&length(size0[-1])>0){
        No<-.RM.EMPTY.NODE(No,size0[-1])
        igraph.index.attr[[i]]<-igraph.index.attr[[i]][-size0[-1],]
        rownames(igraph.index.attr[[i]])<-c(1:nrow(igraph.index.attr[[i]]))
        size0<-size0[1]
      }

      #match cluster.distribution to each vertex
      pie.values<-list()
      for(j in 2:length(igraph.index.attr[[i]]$seq.number)){
        pie.values[[j]]<-cluster.distribution[[igraph.index.attr[[i]]$seq.number[j]]]
      }
      pie.values.list[[i]]<-pie.values


      #write graph list
      g<-igraph::graph(No,directed = T)
      g$layout<-igraph::layout_as_tree
      igraph::V(g)$label<-igraph.index.attr[[i]]$size
      igraph::V(g)$size<-igraph.index.attr[[i]]$adj_size
      igraph::V(g)$label.dist<-3

      igraph::V(g)$shape<-"pie"

      if(length(size0)>1){
        for(k in 1:length(size0)){
          igraph::V(g)$shape[[size0[k]]]<-"circle"
        }
      }
      if(length(size0)==1 ){
        igraph::V(g)$shape[[size0]]<-"circle"
      }

      igraph::V(g)$pie.color<-list(colors)
      igraph::V(g)$color<-"gray"
      igraph::V(g)$color[1]<-"black"
      ###!!!
      igraph::V(g)$pie<-pie.values.list[[i]]

      igraph_list[[i]]<-g

    }

    if (length(igraph.index[[i]])==1){
      igraph.index.attr[[i]]<-size.of.vertex$size[which(size.of.vertex$seq.number==igraph.index[[i]])]#size
      igraph.index.jr[[i]]<-data.frame(igraph.index[[i]])
      pie.values.list[[i]]<-list(cluster.distribution[igraph.index[[i]]][[1]])

      g<-igraph::graph(c(1,1), edges = NULL)
      g$layout<-igraph::layout_as_tree
      igraph::V(g)$size<-igraph.index.attr[[i]]/30+20
      igraph::V(g)$label<-igraph.index.attr[[i]]
      igraph::V(g)$label.dist<-5
      igraph::V(g)$shape<-"pie"
      igraph::V(g)$pie.color<-list(colors)
      igraph::V(g)$pie<-pie.values.list[[i]]

      igraph_list[[i]]<-g
    }
  }
  return(igraph_list)
}

#' Get clone network igraphs without empty mode. Empty node represents the 'extincted' sequences, that are not in any living cell but once existed.
#' @param history The dataframe 'history' from the simulation output.
#' @param igraph.index The list 'igraph.index' from the simulation output.
#' @param empty.node If TRUE, there will be empty node in igraph. if FALSE, the empty node will be deleted.
#' @return  a list of clone network igraphs without empty mode. 
#' @export
no.empty.node<-function(history,igraph.index){
  #require(stats,reshape2)
    igraph.index.attr<-list()
    igraph.index.jr<-list()
    igraph_list<-list()
    igraph_list_trans<-list()
    pie_values_trans_list<-list()
    colors<-c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3" ,"#A6D854")
    pie.values.list<-list()
    pie.values<-list()
    #igraph
    size.of.vertex<-as.data.frame(stats::aggregate(barcode.history~seq.number,history,length,na.action = stats::na.pass))
    #NA count as 1
    size.of.vertex1<-as.data.frame(stats::aggregate(barcode.history~seq.number,history,length))
    size.of.vertex<-merge(size.of.vertex1,size.of.vertex,all=T, by="seq.number")[,-3]
    size.of.vertex[!(size.of.vertex$seq.number %in% size.of.vertex1$seq.number),2]<-0
    #NA count as 0
    colnames(size.of.vertex)<-c("seq.number","size")
    rm(size.of.vertex1)

    #isotype distribution list
    history$isotype[is.na(history$barcode.history)]<-9999
    history$trans_state_his[is.na(history$barcode.history)]<-9999

    isotype.distribution<-reshape2::dcast(history, seq.number~isotype, value.var ="isotype",length)[,-1]
    if(length(which(colnames(isotype.distribution)==9999))!=0){
      isotype.distribution<-isotype.distribution[,1:ncol(isotype.distribution)-1]
    }
    isotype.name<-colnames(isotype.distribution)
    isotype.distribution <-as.list( as.data.frame(t(isotype.distribution)))
    #transcriptome state distribution list
    trans_state_distribution<-reshape2::dcast(history, seq.number~trans_state_his, value.var ="trans_state_his",length)[,-1]
    if(length(which(colnames(trans_state_distribution)==9999))!=0){
      trans_state_distribution<-trans_state_distribution[,1:ncol(trans_state_distribution)-1]
    }
    trans_state_name<-colnames(trans_state_distribution)
    trans_state_distribution <-as.list( as.data.frame(t(trans_state_distribution)))

    #write igraph list and assign attr
    for (i in 1:length(igraph.index)){
      if (length(igraph.index[[i]])>1){
        #igraph.index[[i]]<-c(0,igraph.index[[i]])
        igraph.index.attr[[i]]<-data.frame(igraph.index[[i]])
        colnames(igraph.index.attr[[i]])<-"seq.number"

        igraph.index.jr[[i]]<-data.frame(igraph.index[[i]])
        colnames(igraph.index.jr[[i]])<-"seq.number"

        #delete duplicated for vertex attribute
        igraph.index.attr[[i]]<-subset(igraph.index.attr[[i]],!duplicated(igraph.index.attr[[i]]$seq.number))
        #size of vertex
        igraph.index.attr[[i]]<-merge(size.of.vertex,igraph.index.attr[[i]],by="seq.number",all.y=T)
        igraph.index.attr[[i]]$size[1]<-0
        igraph.index.attr[[i]]$adj_size<-(igraph.index.attr[[i]]$size*60/(range(igraph.index.attr[[i]]$size)[2]-range(igraph.index.attr[[i]]$size)[1]))^0.5+15 #visual size of vertex a number between 20 and 60
        #size of empty
        igraph.index.attr[[i]]$adj_size<-replace(igraph.index.attr[[i]]$adj_size, igraph.index.attr[[i]]$size==0, 10)

        #simplified index
        igraph.index.jr[[i]]$No<-match(x=igraph.index.jr[[i]]$seq.number,table=igraph.index.attr[[i]]$seq.number)

        #find empty
        size0<-which(igraph.index.attr[[i]]$size==0)

        No<-igraph.index.jr[[i]]$No
        #option: no empty node

        No<-.RM.EMPTY.NODE(No,size0[-1])
        if(length(size0[-1])>0){
          igraph.index.attr[[i]]<-igraph.index.attr[[i]][-size0[-1],]
        }

        rownames(igraph.index.attr[[i]])<-c(1:nrow(igraph.index.attr[[i]]))

        #match isotype.distribution to each vertex
        pie.values<-list()
        for(j in 2:length(igraph.index.attr[[i]]$seq.number)){
          pie.values[[j]]<-isotype.distribution[[igraph.index.attr[[i]]$seq.number[j]]]
        }
        pie.values.list[[i]]<-pie.values

        #match trans state distribution to each vertex
        pie_values_trans<-list()
        for(j in 2:length(igraph.index.attr[[i]]$seq.number)){
          pie_values_trans[[j]]<-trans_state_distribution[[igraph.index.attr[[i]]$seq.number[j]]]
        }
        pie_values_trans_list[[i]]<-pie_values_trans

        #write graph list
        g<-igraph::graph(No,directed = T)
        g$layout<-igraph::layout_as_tree
        igraph::V(g)$label<-igraph.index.attr[[i]]$size
        igraph::V(g)$size<-igraph.index.attr[[i]]$adj_size
        igraph::V(g)$label.dist<-3

        igraph::V(g)$shape<-"pie"
        igraph::V(g)$shape[1]<-"circle"
        igraph::V(g)$pie.color<-list(colors)
        igraph::V(g)$color[1]<-"black"
        p<-g
        ###!!!
        igraph::V(g)$pie<-pie.values.list[[i]]
        igraph::V(p)$pie<-pie_values_trans_list[[i]]

        igraph_list[[i]]<-g
        igraph_list_trans[[i]]<-p
      }

      if (length(igraph.index[[i]])==1){
        igraph.index.attr[[i]]<-size.of.vertex$size[which(size.of.vertex$seq.number==igraph.index[[i]])]#size
        igraph.index.jr[[i]]<-data.frame(igraph.index[[i]])
        pie.values.list[[i]]<-list(isotype.distribution[igraph.index[[i]]][[1]])
        pie_values_trans_list[[i]]<-list(trans_state_distribution[igraph.index[[i]]][[1]])

        g<-igraph::graph(c(1,1), edges = NULL)
        g$layout<-igraph::layout_as_tree
        igraph::V(g)$size<-igraph.index.attr[[i]]/30+20
        igraph::V(g)$label<-igraph.index.attr[[i]]
        igraph::V(g)$label.dist<-5
        igraph::V(g)$shape<-"pie"
        igraph::V(g)$pie.color<-list(colors)
        p<-g
        igraph::V(g)$pie<-pie.values.list[[i]]
        igraph::V(p)$pie<-pie_values_trans_list[[i]]

        igraph_list[[i]]<-g
        igraph_list_trans[[i]]<-p
      }
    }
    output.list<-list(igraph_list,igraph_list_trans)
    names(output.list)<-c("igraph_list_iso","igraph_list_trans")
    return(output.list)
}

