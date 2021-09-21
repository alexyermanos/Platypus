## ----setup, include=FALSE-----------------------------------------------------
#to avoid of the dplyr summarise warnings
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(stats, warn.conflicts = F)

knitr::opts_chunk$set(fig.width=7, fig.height=7) 

## ---- eval = FALSE, include=FALSE---------------------------------------------
#  #this is for sourcing all functions under development and is not included in the knitting
#  source_code_dir <- "C:/Dokumente usw/Master/Reddy Lab/1_thesis/3_code/4_DB project/PlatypusDB_admin/R"
#  file_path_vec <- list.files(source_code_dir, full.names = T)
#  for(f_path in file_path_vec){
#    print(f_path)
#    tryCatch({source(f_path)}, error = function(e){e})
#  }
#  
#  knitr::opts_chunk$set(fig.width=7, fig.height=7)
#  
#  gc()
#  
#  load("C:/Dokumente usw/Master/Reddy Lab/1_thesis/3_code/4_DB project/PlatypusDB_dev_VK/R/ABForests_sysdata.rda")
#  
#  load("C:/Dokumente usw/Master/Reddy Lab/1_thesis/3_code/4_DB project/PlatypusDB_dev_VK/data/Bcell_sequences_example_tree.rda")
#  file1 <- Bcell_sequences_example_tree
#  
#  load("C:/Dokumente usw/Master/Reddy Lab/1_thesis/3_code/4_DB project/PlatypusDB_dev_VK/data/Bcell_tree_2.rda")
#  file2 <- Bcell_tree_2

## ---- eval = FALSE ,  warning=FALSE,message=FALSE-----------------------------
#  library(ggplot2)
#  library(scales)
#  library(gridExtra)
#  library(grid)

## ---- eval = FALSE ,  include=TRUE, fig.align="center", fig.cap=c("Flowchart summarizing the workflow of AbNetForests. The link with the R package Platypus (Yermanos et al. 2021) is also depicted."), echo=FALSE,out.width='0.75\\linewidth'----
#  knitr::include_graphics("general_sketch.pdf")

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  # Program parameters
#  DAG<-TRUE
#  clonal_frequency<-TRUE
#  scaleByClonalFreq<-TRUE
#  weight<-TRUE
#  tie_flag<-'rand'
#  scaleBybetweenness<-FALSE
#  scaleByclocloseness_metr<-FALSE
#  topdf<-TRUE
#  custom_mat<-NULL
#  random.seed<-3
#  alg_opt<-'naive'
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  #"new" is an example dataset loaded with the package (See ABForests_sysdata)
#  
#  graphs<-AbForests_AntibodyForest(new,FALSE,NULL,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,"isotype",random.seed,alg_opt,NULL)
#  

## ---- eval = FALSE , echo=T, results='hide'-----------------------------------
#  
#  AbForests_PlotGraphs(graphs,NULL,topdf,"Networks")
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs[14][[1]]$network,vertex.label.cex=1+(0.01*graphs[14][[1]]$cells.per.variant),vertex.size=12+(1.2*graphs[14][[1]]$cells.per.variant))
#  
#  legend(graphs[14][[1]]$legend[[1]],bty=graphs[14][[1]]$legend[[2]],legend=graphs[14][[1]]$legend[[3]],fill=graphs[14][[1]]$legend[[4]],border=graphs[14][[1]]$legend[[5]],inset=graphs[14][[1]]$legend[[6]],xpd=graphs[14][[1]]$legend[[7]],title=graphs[14][[1]]$legend$title,cex=graphs[14][[1]]$legend[[10]]*0.7)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  metrics<-AbForests_ForestMetrics(graphs,DAG,clonal_frequency,scaleByClonalFreq,weight,tie_flag,"isotype")
#  

## ---- eval = FALSE , echo=T, results='hide'-----------------------------------
#  
#  AbForests_PlotGraphs(metrics,6,topdf,"Metrics")
#  

## ---- eval = FALSE,  fig.show='hold',echo=T, results='hide'-------------------
#  
#  graphs_cluster<-AbForests_AntibodyForest(new,FALSE,NULL,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,"cluster",random.seed,alg_opt,NULL)
#  

## ---- eval = FALSE , echo=T, results='hide'-----------------------------------
#  
#  AbForests_PlotGraphs(graphs_cluster,NULL,topdf,"Networks_cluster")
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_cluster[14][[1]]$network,vertex.label.cex=1+(0.01*graphs_cluster[14][[1]]$cells.per.variant),vertex.size=12+(1.2*graphs_cluster[14][[1]]$cells.per.variant))
#  
#  legend(graphs_cluster[14][[1]]$legend[[1]],bty=graphs_cluster[14][[1]]$legend[[2]],legend=graphs_cluster[14][[1]]$legend[[3]],fill=graphs_cluster[14][[1]]$legend[[4]],border=graphs_cluster[14][[1]]$legend[[5]],inset=graphs_cluster[14][[1]]$legend[[6]],xpd=graphs_cluster[14][[1]]$legend[[7]],title=graphs_cluster[14][[1]]$legend$title,cex=graphs_cluster[14][[1]]$legend[[10]]*0.35)
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  head(new[[1]],3)

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  file1<-data("Bcell_sequences_example_tree", package = "Platypus")
#  file1<-get(file1)
#  head(file1,3)

## ---- eval = FALSE,  fig.show='hold'------------------------------------------
#  file2<-data(Bcell_tree_2, package = "Platypus")
#  file2<-get(file2)
#  head(file2,3)

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  files<-c()
#  files<-append(files,file1)
#  files<-append(files,file2)

## ---- eval = FALSE------------------------------------------------------------
#  graphs_csv<-AbForests_AntibodyForest(files,TRUE,NULL,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,"isotype",random.seed,alg_opt,0)

## ---- eval = FALSE------------------------------------------------------------
#  igraph::plot.igraph(graphs_csv[1][[1]]$network,vertex.label.cex=1+(0.01*graphs_csv[1][[1]]$cells.per.variant),vertex.size=8+(1.2*graphs_csv[1][[1]]$cells.per.variant))
#  
#  legend(graphs_csv[1][[1]]$legend[[1]],bty=graphs_csv[1][[1]]$legend[[2]],legend=graphs_csv[1][[1]]$legend[[3]],fill=graphs_csv[1][[1]]$legend[[4]],border=graphs_csv[1][[1]]$legend[[5]],inset=graphs_csv[1][[1]]$legend[[6]],xpd=graphs_csv[1][[1]]$legend[[7]],title=graphs_csv[1][[1]]$legend$title,cex=graphs_csv[1][[1]]$legend[[10]]*0.7)

## ---- eval = FALSE------------------------------------------------------------
#  
#  graphs_csv<-AbForests_AntibodyForest(files,TRUE,NULL,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,"isotype",random.seed,alg_opt,1)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_csv[1][[1]]$network,vertex.label.cex=1+(0.01*graphs_csv[1][[1]]$cells.per.variant),vertex.size=8+(1.2*graphs_csv[1][[1]]$cells.per.variant))
#  
#  legend(graphs_csv[1][[1]]$legend[[1]],bty=graphs_csv[1][[1]]$legend[[2]],legend=graphs_csv[1][[1]]$legend[[3]],fill=graphs_csv[1][[1]]$legend[[4]],border=graphs_csv[1][[1]]$legend[[5]],inset=graphs_csv[1][[1]]$legend[[6]],xpd=graphs_csv[1][[1]]$legend[[7]],title=graphs_csv[1][[1]]$legend$title,cex=graphs_csv[1][[1]]$legend[[10]]*0.7)
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  total_lineages<-length(new)
#  total_lineages
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  df_total_lineages <- data.frame(total_lineages=total_lineages,colors=c('black'))
#  df_total_lineages$x_values <-c(rep("1",each=nrow(df_total_lineages)))
#  df_total_lineages$colors<-as.factor(df_total_lineages$colors)
#  
#  p_total_lineages<-ggplot(df_total_lineages, aes(x = x_values, y =total_lineages ,fill=colors)) + ylab("Number of clonal lineages")+
#   coord_cartesian() +
#    scale_x_discrete(expand = expansion(add=c(0.6,0.5)))+
#    scale_y_continuous(expand = c(0,0)) +
#    geom_bar(width=0.5,stat="identity", position="dodge",colour = "black",show.legend = FALSE)+
#    theme(aspect.ratio =1,
#          axis.title.x=element_blank(),
#          axis.text.x=element_blank(),
#          axis.ticks.x=element_blank(),
#          axis.ticks.length=unit(.3, "cm"),
#          axis.title.y = element_text(size = rel(1), angle = 90),
#          axis.text.y = element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
#          axis.text = element_text(size = 8),
#          panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          panel.background = element_blank(),
#          axis.line = element_line(colour = "black"))+
#          scale_fill_manual(values = levels(df_total_lineages$colors))
#  
#  p_total_lineages
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  bcells_per_clonal_lineage<-unname(unlist(sapply(new,function(z) nrow(z))))
#  df<-data.frame(Combinations=1:total_lineages, Freq=bcells_per_clonal_lineage)
#  p_cells_all<-ggplot(df,aes_string(x = 'Combinations', y = 'Freq')) + xlab("Clonal lineages") + ylab("Number of cells")+geom_bar(stat="identity", position="dodge",colour = "black",show.legend = FALSE)+
#    coord_cartesian(expand = FALSE )+
#    theme(aspect.ratio=1,axis.title.y = element_text(size = rel(1), angle = 90),
#          axis.title.x = element_text(size = rel(1), angle = 0),
#          axis.line = element_line(colour = "black"),
#          panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          axis.ticks.length=unit(.3, "cm"),
#          axis.text = element_text(size = 8),
#          axis.text.x = element_text(angle = 0, vjust = 0.5,colour = "black",face="bold"),
#          axis.ticks.x = element_blank(),
#          axis.text.y = element_text(angle = 0, vjust = 0.5,colour = "black",face="bold"),
#          panel.background = element_blank(),legend.title=element_text(size=10,face="bold")) + scale_fill_discrete(name = "Vertex Degrees")
#  
#  p_cells_all
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  mean_bcells_per_clonal_lineage<-round(mean(bcells_per_clonal_lineage))
#  mean_bcells_per_clonal_lineage
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  std_bcells_per_clonal_lineage<-round(sd(bcells_per_clonal_lineage))
#  std_bcells_per_clonal_lineage
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  unique_Ab_variants_per_clonal_lineage<- AbForests_UniqueAntibodyVariants(new,"isotype",custom_mat,tie_flag,weight,random.seed,alg_opt,NULL)
#  df_t<-data.frame(Combinations=1:total_lineages, Freq=unique_Ab_variants_per_clonal_lineage)
#  p_variants_all<-ggplot(df_t,aes_string(x = 'Combinations', y = 'Freq')) + xlab("Clonal lineages") + ylab("Number of unique \n antibody variants")+geom_bar(stat="identity", position="dodge",colour = "black",show.legend = FALSE)+
#    coord_cartesian(expand = FALSE )+
#    theme(aspect.ratio=1,axis.title.y = element_text(size = rel(1), angle = 90),
#          axis.title.x = element_text(size = rel(1), angle = 0),
#          axis.text = element_text(size = 8),axis.line = element_line(colour = "black"),
#          panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          axis.text.x = element_text(angle = 0, vjust = 0.5,colour = "black",face="bold"),
#          axis.ticks.x = element_blank(),
#          axis.text.y = element_text(angle = 00, vjust = 0.5,colour = "black",face="bold"),
#          panel.background = element_blank(),legend.title=element_text(size=10,face="bold")) + scale_fill_discrete(name = "Vertex Degrees")
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  df_mean_cells_variants <- data.frame(Freq=c(bcells_per_clonal_lineage,unique_Ab_variants_per_clonal_lineage))
#  df_mean_cells_variants$sample_name <-c(rep("Cells",each=length(bcells_per_clonal_lineage)),rep("Variants",each=length(unique_Ab_variants_per_clonal_lineage)))
#  df_mean_cells_variants$colors <-as.vector(c(rep("black",each=length(bcells_per_clonal_lineage)),rep("gray85",each=length(unique_Ab_variants_per_clonal_lineage))))
#  df_mean_cells_variants$colors<-as.factor(df_mean_cells_variants$colors)
#  
#  p_mean_cells_variants<-ggplot(df_mean_cells_variants, aes(sample_name, Freq ,fill=colors)) + ylab("Number per clonal lineage")+ coord_cartesian() +
#    geom_errorbar(stat = "summary", fun.data = "mean_se",
#                  fun.args = list(mult = 1),
#                  position ="dodge", width = 0.2) +
#    geom_bar(stat = "summary", fun = "mean",
#             position = position_dodge(width = 0.2), show.legend = FALSE) +
#    scale_y_continuous(expand = c(0,0)) +
#    scale_x_discrete(name = "sample_name") +
#    theme(aspect.ratio = 1,axis.title.x=element_blank(),
#      axis.title.y = element_text(size = rel(1), angle = 90),
#      axis.text = element_text(size = 8),axis.line = element_line(colour = "black"),
#          panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          axis.text.x = element_text(angle = 0, vjust = 0.9,colour = "black",face="bold"),
#          axis.ticks.x = element_blank(),
#      axis.text.y = element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
#          panel.background = element_blank(),legend.title=element_text(size=10,face="bold"))+scale_fill_manual(values = levels(df_mean_cells_variants$colors))
#  
#  p_mean_cells_variants

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  mean_unique_Ab_variants_per_clonal_lineage<-round(mean(unique_Ab_variants_per_clonal_lineage))
#  mean_unique_Ab_variants_per_clonal_lineage
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  std_unique_Ab_variants_per_clonal_lineage<-round(sd(unique_Ab_variants_per_clonal_lineage))
#  std_unique_Ab_variants_per_clonal_lineage
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  pos<-which(bcells_per_clonal_lineage == max(bcells_per_clonal_lineage[bcells_per_clonal_lineage>0]))
#  head(new[[pos]],3)
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  pos<-which(unique_Ab_variants_per_clonal_lineage == max(unique_Ab_variants_per_clonal_lineage[unique_Ab_variants_per_clonal_lineage>0]))
#  head(new[[pos]],3)
#  

## ---- eval = FALSE, warning=FALSE---------------------------------------------
#  
#  p_total_lineages <- arrangeGrob(p_total_lineages, top = textGrob("A", x = unit(0, "npc"),y= unit(1, "npc"), just=c("left","top"),                                 gp=gpar(col="black", fontsize=12, fontfamily="Times")))
#  
#  p_cells_all <- arrangeGrob(p_cells_all, top = textGrob("B", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"),
#  gp=gpar(col="black", fontsize=12, fontfamily="Times")))
#  
#  p_variants_all <- arrangeGrob(p_variants_all, top = textGrob("C", x = unit(0, "npc"),y  = unit(1, "npc"), just=c("left","top"),
#  gp=gpar(col="black", fontsize=12, fontfamily="Times")))
#  
#  p_mean_cells_variants <- arrangeGrob(p_mean_cells_variants, top = textGrob("D", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black",fontsize=12, fontfamily="Times")))
#  
#  grid.arrange(p_total_lineages,p_cells_all,
#              p_variants_all,p_mean_cells_variants,
#              ncol=2)
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  repert<- AbForests_SubRepertoiresByCells(new)
#  

## ---- eval = FALSE ,  fig.show='hold'-----------------------------------------
#  
#  total_lineages_per_isotype<-unlist(lapply(repert, function(x) length(x)))
#  total_lineages_per_isotype
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  df_clon_lineages_isotype <- data.frame(Strategies= c("IgG","IgA","IgM"),Freq=unname(total_lineages_per_isotype[1:3]),colors=c("green1","red","purple"))
#  df_clon_lineages_isotype$colors<-as.factor(df_clon_lineages_isotype$colors)
#  
#  p_clon_lineages_isotype<-ggplot(df_clon_lineages_isotype, aes_string(x = 'Strategies', y = 'Freq', fill='colors')) + ylab("# of clonal lineages \n before filetering")+
#    geom_bar(width=0.6,stat="identity", position="dodge",colour = "black",show.legend = FALSE)+
#    coord_cartesian()+
#    scale_x_discrete(expand = expansion(add=c(0.6,0.5)))+
#    scale_y_continuous(expand = c(0,0)) +
#    theme(aspect.ratio = 1,axis.title.x=element_blank(),
#          axis.title.y = element_text(size = rel(1), angle = 90),
#          axis.text = element_text(size = 8),axis.line = element_line(colour = "black"),
#          panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          axis.text.x = element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
#          axis.text.y = element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
#          panel.background = element_blank(),legend.title=element_text(size=8,face="bold"))+scale_fill_manual(values = levels(df_clon_lineages_isotype$colors))
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  list6<-AbForests_ConvertStructure(repert,"isotype",NULL)
#  rev_rep<-AbForests_RemoveNets(list6,"isotype",custom_mat,tie_flag,weight,NULL,6,random.seed,alg_opt)
#  remained_lineages_per_isotype<-unlist(lapply(rev_rep, function(x) length(x)))
#  remained_lineages_per_isotype
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  df_clon_lineages_isotype_after <- data.frame(Strategies=c("IgG","IgA","IgM"), Freq=unname(remained_lineages_per_isotype[1:3]),colors=c("green1","red","purple"))
#  df_clon_lineages_isotype_after$colors <- as.factor(df_clon_lineages_isotype_after$colors)
#  
#  p_clon_lineages_isotype_after<-ggplot(df_clon_lineages_isotype_after, aes_string(x = 'Strategies', y = 'Freq',fill='colors')) + ylab("# of clonal lineages \n after filtering")+
#    geom_bar(width=0.6,stat="identity", position="dodge",colour = "black",show.legend = FALSE)+
#    coord_cartesian()+
#    scale_x_discrete(expand = expansion(add=c(0.6,0.5)))+
#    scale_y_continuous(expand = c(0,0)) +
#    theme(aspect.ratio=1, axis.title.x=element_blank(),
#          axis.title.y = element_text(size = rel(1), angle = 90),
#          axis.text = element_text(size = 8),axis.line = element_line(colour = "black"),
#          panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          axis.text.x = element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
#          axis.text.y = element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
#          panel.background = element_blank(),legend.title=element_text(size=8,face="bold"))+scale_fill_manual(values = levels(df_clon_lineages_isotype_after$colors))
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  percent_remained<-unlist(mapply("/",remained_lineages_per_isotype[1:3],total_lineages_per_isotype[1:3],SIMPLIFY = FALSE))
#  
#  df_percent_remained <- data.frame(Strategies=c("IgG","IgA","IgM"),Freq=unlist(unname(percent_remained)),colors=c("green1","red","purple"))
#  
#  df_percent_remained[] <- lapply(df_percent_remained, function(x) if(is.numeric(x)) x*100 else x)
#  df_percent_remained$colors <- factor(df_percent_remained$colors, levels = df_percent_remained$colors)
#  p_rem<-ggplot(df_percent_remained, aes(x = Strategies, y = Freq, fill= colors)) + ylab("% of clonal lineages \n after filtering")+
#    geom_bar(width=0.5,stat="identity", position="dodge",colour = "black",show.legend = FALSE)+
#    coord_cartesian()+
#    scale_x_discrete(expand = expansion(add=c(0.6,0.5)))+
#    scale_y_continuous(expand = c(0,0)) +
#    theme(aspect.ratio=1,axis.line = element_line(colour = "black"),
#          axis.title.x=element_blank(),
#          panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          axis.title.y = element_text(size = rel(1), angle = 90),
#          axis.text.y = element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
#          axis.text = element_text(size = 8),
#          axis.text.x = element_text(angle = 0, vjust = 0.5,colour = "black",face="bold"),
#          panel.background = element_blank(),legend.title=element_text(size=8,face="bold")) + scale_fill_manual(values = levels(df_percent_remained$colors))
#  

## ---- eval = FALSE, warning=FALSE---------------------------------------------
#  
#  p_clon_lineages_isotype <- arrangeGrob(p_clon_lineages_isotype, top = textGrob("A", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=12, fontfamily="Times")))
#  
#  p_clon_lineages_isotype_after <- arrangeGrob(p_clon_lineages_isotype_after, top = textGrob("B", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"), gp=gpar(col="black",fontsize=12, fontfamily="Times")))
#  
#  p_rem <- arrangeGrob(p_rem, top = textGrob("C", x = unit(0, "npc")
#  , y = unit(1, "npc"), just=c("left","top"),  gp=gpar(col="black",    fontsize=12, fontfamily="Times")))
#  
#  grid.arrange(p_clon_lineages_isotype,
#              p_clon_lineages_isotype_after,p_rem, ncol=2)
#  

## ---- eval = FALSE------------------------------------------------------------
#  graphs_IGHG<-AbForests_AntibodyForest(NULL,FALSE,rev_rep$list_IGHG,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,"isotype",random.seed,alg_opt,NULL)

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_IGHG[11][[1]]$network,vertex.label.cex=1+(0.0001*graphs_IGHG[11][[1]]$cells.per.variant),vertex.size=17+(0.8*igraph::V(graphs_IGHG[11][[1]]$network)$size),layout=0.05*graphs_IGHG[11][[1]]$network$layout)
#  legend(graphs_IGHG[11][[1]]$legend[[1]],bty = graphs_IGHG[11][[1]]$legend[[2]],legend=graphs_IGHG[11][[1]]$legend[[3]],fill=graphs_IGHG[11][[1]]$legend[[4]],border=graphs_IGHG[11][[1]]$legend[[5]],inset=graphs_IGHG[11][[1]]$legend[[6]],xpd=graphs_IGHG[11][[1]]$legend[[7]],title=graphs_IGHG[11][[1]]$legend$title,cex=graphs_IGHG[11][[1]]$legend[[10]]*0.7)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_IGHG[104][[1]]$network,vertex.label.cex=1+(0.0001*graphs_IGHG[104][[1]]$cells.per.variant),vertex.size=17+(0.8*igraph::V(graphs_IGHG[104][[1]]$network)$size),layout=0.05*graphs_IGHG[104][[1]]$network$layout)
#  legend(graphs_IGHG[104][[1]]$legend[[1]],bty = graphs_IGHG[104][[1]]$legend[[2]],legend=graphs_IGHG[104][[1]]$legend[[3]],fill=graphs_IGHG[104][[1]]$legend[[4]],border=graphs_IGHG[104][[1]]$legend[[5]],inset=graphs_IGHG[104][[1]]$legend[[6]],xpd=graphs_IGHG[104][[1]]$legend[[7]],title=graphs_IGHG[104][[1]]$legend$title,cex=graphs_IGHG[104][[1]]$legend[[10]]*0.7)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  graphs_IGHM<-AbForests_AntibodyForest(NULL,FALSE,rev_rep$list_IGHM,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,"isotype",random.seed,alg_opt,NULL)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_IGHM[1][[1]]$network,vertex.label.cex=1+(0.0001*graphs_IGHM[1][[1]]$cells.per.variant),vertex.size=17+(0.8*igraph::V(graphs_IGHM[1][[1]]$network)$size),layout=0.05*graphs_IGHM[1][[1]]$network$layout)
#  legend(graphs_IGHM[1][[1]]$legend[[1]],bty = graphs_IGHM[1][[1]]$legend[[2]],legend=graphs_IGHM[1][[1]]$legend[[3]],fill=graphs_IGHM[1][[1]]$legend[[4]],border=graphs_IGHM[1][[1]]$legend[[5]],inset=graphs_IGHM[1][[1]]$legend[[6]],xpd=graphs_IGHM[1][[1]]$legend[[7]],title=graphs_IGHM[1][[1]]$legend$title,cex=graphs_IGHM[1][[1]]$legend[[10]]*0.7)
#  

## ---- eval = FALSE------------------------------------------------------------
#  #Scaling of nodes by vertex betweenness and labeling based on traversal order of nodes in the naive algorithm
#  graphs_IGHM_new<-AbForests_AntibodyForest(NULL,FALSE,rev_rep$list_IGHM,custom_mat,clonal_frequency=FALSE,scaleByClonalFreq=FALSE,weight,tie_flag,scaleBybetweenness=TRUE,scaleByclocloseness_metr=FALSE,"isotype",random.seed,alg_opt,NULL)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_IGHM_new[1][[1]]$network,vertex.label.cex=1+(0.0001*graphs_IGHM_new[1][[1]]$cells.per.variant),vertex.size=17+(0.8*igraph::V(graphs_IGHM_new[1][[1]]$network)$size),layout=0.05*graphs_IGHM_new[1][[1]]$network$layout)
#  legend(graphs_IGHM_new[1][[1]]$legend[[1]],bty = graphs_IGHM_new[1][[1]]$legend[[2]],legend=graphs_IGHM_new[1][[1]]$legend[[3]],fill=graphs_IGHM_new[1][[1]]$legend[[4]],border=graphs_IGHM_new[1][[1]]$legend[[5]],inset=graphs_IGHM_new[1][[1]]$legend[[6]],xpd=graphs_IGHM_new[1][[1]]$legend[[7]],title=graphs_IGHM_new[1][[1]]$legend$title,cex=graphs_IGHM_new[1][[1]]$legend[[10]]*0.7)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  graphs_IGHG_2<-AbForests_AntibodyForest(NULL,FALSE,rev_rep$list_IGHG,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,"isotype",random.seed,alg_opt="two-step",NULL)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_IGHG_2[104][[1]]$network,vertex.label.cex=1+(0.0001*graphs_IGHG_2[104][[1]]$cells.per.variant),vertex.size=17+(0.8*igraph::V(graphs_IGHG_2[104][[1]]$network)$size),layout=0.05*graphs_IGHG_2[104][[1]]$network$layout)
#  legend(graphs_IGHG_2[104][[1]]$legend[[1]],bty = graphs_IGHG_2[104][[1]]$legend[[2]],legend=graphs_IGHG_2[104][[1]]$legend[[3]],fill=graphs_IGHG_2[104][[1]]$legend[[4]],border=graphs_IGHG_2[104][[1]]$legend[[5]],inset=graphs_IGHG_2[104][[1]]$legend[[6]],xpd=graphs_IGHG_2[104][[1]]$legend[[7]],title=graphs_IGHG_2[104][[1]]$legend$title,cex=graphs_IGHG_2[104][[1]]$legend[[10]]*0.7)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  graphs_IGHG_close_germ<-AbForests_AntibodyForest(NULL,FALSE,rev_rep$list_IGHG,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag="close_to_germ",scaleBybetweenness,scaleByclocloseness_metr,"isotype",random.seed,alg_opt="two-step",NULL)
#  

## ---- eval = FALSE------------------------------------------------------------
#  igraph::plot.igraph(graphs_IGHG_close_germ[104][[1]]$network,vertex.label.cex=1+(0.0001*graphs_IGHG_close_germ[104][[1]]$cells.per.variant),vertex.size=17+(0.8*igraph::V(graphs_IGHG_close_germ[104][[1]]$network)$size),layout=0.05*graphs_IGHG_close_germ[104][[1]]$network$layout)
#  legend(graphs_IGHG_close_germ[104][[1]]$legend[[1]],bty = graphs_IGHG_close_germ[104][[1]]$legend[[2]],legend=graphs_IGHG_close_germ[104][[1]]$legend[[3]],fill=graphs_IGHG_close_germ[104][[1]]$legend[[4]],border=graphs_IGHG_close_germ[104][[1]]$legend[[5]],inset=graphs_IGHG_close_germ[104][[1]]$legend[[6]],xpd=graphs_IGHG_close_germ[104][[1]]$legend[[7]],title=graphs_IGHG_close_germ[104][[1]]$legend$title,cex=graphs_IGHG_close_germ[104][[1]]$legend[[10]]*0.7)

## ---- eval = FALSE------------------------------------------------------------
#  
#  graphs_IGHG_full<-AbForests_AntibodyForest(NULL,FALSE,rev_rep$list_IGHG,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag="full",scaleBybetweenness,scaleByclocloseness_metr,"isotype",random.seed,alg_opt="two-step",NULL)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_IGHG_full[104][[1]]$network,vertex.label.cex=1+(0.0001*graphs_IGHG_full[104][[1]]$cells.per.variant),vertex.size=17+(0.8*igraph::V(graphs_IGHG_full[104][[1]]$network)$size),layout=0.05*graphs_IGHG_full[104][[1]]$network$layout)
#  legend(graphs_IGHG_full[104][[1]]$legend[[1]],bty = graphs_IGHG_full[104][[1]]$legend[[2]],legend=graphs_IGHG_full[104][[1]]$legend[[3]],fill=graphs_IGHG_full[104][[1]]$legend[[4]],border=graphs_IGHG_full[104][[1]]$legend[[5]],inset=graphs_IGHG_full[104][[1]]$legend[[6]],xpd=graphs_IGHG_full[104][[1]]$legend[[7]],title=graphs_IGHG_full[104][[1]]$legend$title,cex=graphs_IGHG_full[104][[1]]$legend[[10]]*0.7)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  metrics_igg<-AbForests_ForestMetrics(graphs_IGHG,DAG,clonal_frequency,scaleByClonalFreq,weight,tie_flag,"isotype")
#  
#  metrics_igg[[1]]$`Mean Total path length from germline`
#  

## ---- eval = FALSE------------------------------------------------------------
#  metrics_igg[[1]]$`Plot: Distribution of isotypes`

## ---- eval = FALSE , echo=T, results='hide'-----------------------------------
#  output<-AbForests_CompareForests(graphs_IGHG,graphs_IGHM,DAG,clonal_frequency,scaleByClonalFreq,weight,tie_flag,"isotype")

## ---- eval = FALSE------------------------------------------------------------
#  
#  graphs_custom<-AbForests_AntibodyForest(new,FALSE,NULL,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,"VDJ_jgene",random.seed,alg_opt,NULL)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_custom[14][[1]]$network,vertex.label.cex=1+(0.0001*graphs_custom[14][[1]]$cells.per.variant),vertex.size=17+(0.8*igraph::V(graphs_custom[14][[1]]$network)$size),layout=0.05*graphs_custom[14][[1]]$network$layout)
#  legend(graphs_custom[14][[1]]$legend[[1]],bty = graphs_custom[14][[1]]$legend[[2]],legend=graphs_custom[14][[1]]$legend[[3]],fill=graphs_custom[14][[1]]$legend[[4]],border=graphs_custom[14][[1]]$legend[[5]],inset=graphs_custom[14][[1]]$legend[[6]],xpd=graphs_custom[14][[1]]$legend[[7]],title=graphs_custom[14][[1]]$legend$title,cex=graphs_custom[14][[1]]$legend[[10]]*0.7)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  graphs_csv<-AbForests_AntibodyForest(files,TRUE,NULL,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,"LCvgene",random.seed,alg_opt,0)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_csv[1][[1]]$network,vertex.label.cex=1+(0.0001*graphs_csv[1][[1]]$cells.per.variant),vertex.size=17+(0.8*igraph::V(graphs_csv[1][[1]]$network)$size),layout=0.05*graphs_csv[1][[1]]$network$layout)
#  legend(graphs_csv[1][[1]]$legend[[1]],bty = graphs_csv[1][[1]]$legend[[2]],legend=graphs_csv[1][[1]]$legend[[3]],fill=graphs_csv[1][[1]]$legend[[4]],border=graphs_csv[1][[1]]$legend[[5]],inset=graphs_csv[1][[1]]$legend[[6]],xpd=graphs_csv[1][[1]]$legend[[7]],title=graphs_csv[1][[1]]$legend$title,cex=graphs_csv[1][[1]]$legend[[10]]*0.7)
#  

## ---- eval = FALSE------------------------------------------------------------
#  repert<-AbForests_SubRepertoiresByUniqueSeq(new,"isotype",custom_mat,tie_flag,weight,random.seed,alg_opt,0)

## ---- eval = FALSE------------------------------------------------------------
#  
#  list_scv<-AbForests_CsvToDf(files)
#  list6<-AbForests_ConvertStructure(list_scv,"cluster",0)
#  rev_rep<-AbForests_RemoveNets(list6,"cluster",custom_mat,tie_flag,weight,4,NULL,random.seed,alg_opt)
#  graphs_csv_less_4<-AbForests_AntibodyForest(NULL,FALSE,rev_rep,custom_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,"cluster",random.seed,alg_opt,0)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  igraph::plot.igraph(graphs_csv_less_4[1][[1]]$network,vertex.label.cex=1+(0.0001*graphs_csv_less_4[1][[1]]$cells.per.variant),vertex.size=17+(0.8*igraph::V(graphs_csv_less_4[1][[1]]$network)$size),layout=0.05*graphs_csv_less_4[1][[1]]$network$layout)
#  legend(graphs_csv_less_4[1][[1]]$legend[[1]],bty = graphs_csv_less_4[1][[1]]$legend[[2]],legend=graphs_csv_less_4[1][[1]]$legend[[3]],fill=graphs_csv_less_4[1][[1]]$legend[[4]],border=graphs_csv_less_4[1][[1]]$legend[[5]],inset=graphs_csv_less_4[1][[1]]$legend[[6]],xpd=graphs_csv_less_4[1][[1]]$legend[[7]],title=graphs_csv_less_4[1][[1]]$legend$title,cex=graphs_csv_less_4[1][[1]]$legend[[10]]*0.26)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  AbForests_PlyloToMatrix("example.tree")
#  

