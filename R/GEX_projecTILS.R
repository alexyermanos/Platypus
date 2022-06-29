#' Projection of scRNA-seq data into reference single-cell atlas, enabling their celltype annotation based on the single-cell atlas.
#' @param ref_path Path to reference TIL atlas file (ex: c:/Users/.../ref_TILAtlas_mouse_v1.rds). The atlas can be downloaded from the GitHub of ProjecTILs.
#' @param GEX GEX output of the VDJ_GEX_matrix function (VDJ_GEX_matrix[[2]])).
#' @param split_by Optional character vector to specify how the GEX should be split for analysis. This parameter can refer to any column in the GEX. If none is given by the user the analysis will take the whole GEX.
#' @param filtering Logical, if TRUE a filtering is apply which eliminates unwanted cells. By default it is set to FALSE.
#' @param NA_cells Logical, if TRUE the cells not assigned by projecTILs are kept in the bar plot, if FALSE, not assigned cells are filtered out. By default it is set to TRUE.
#' @return Return a list. Element[[1]] is the GEX data frame containing two new columns containing ProjecTILs cell type assignment. Element[[2]] is the output of make.projection function from projecTILs based on the given GEX. Element[[3]] contains a UMAP plot per each groups based on projecTILs assignment. Element[[4]] plots of the fraction of cells with predicted state per cluster. 
#' @export
#' @examples
#' \dontrun{
#' #Without splitting argument, considering the whole VGM
#' output_projecTILS_wohle_VGM<-GEX_projecTILS(ref_path = "c:/Users/.../ref_TILAtlas_mouse_v1.rds", GEX = VGM$GEX,
#' filtering =TRUE)
#' output_projecTILS_wohle_VGM[[3]] #Umap
#' output_projecTILS_wohle_VGM[[4]] #Barplots
#' 
#' #With splitting argument by groups_id
#' output_projecTILS_split_by_group<-GEX_projecTILS(ref_path = "c:/Users/.../ref_TILAtlas_mouse_v1.rds", GEX = VGM$GEX, 
#' filtering =  TRUE, split_by = "group_id", NA_cells = FALSE)
#' output_projecTILS_split_by_group[[3]] #Umap
#' output_projecTILS_split_by_group[[4]] #Barplots
#'}

GEX_projecTILS<-function(ref_path, GEX, split_by, filtering=c(TRUE,FALSE),NA_cells=c(TRUE,FALSE)){
  if(missing(ref_path))stop("Please provide reference path to reference TIL atlas file")
  
  platypus.version <- "It doesn't matter"
  
  if(missing(NA_cells)){
    NA_cells = TRUE
  }
  if(missing(filtering)){filtering = FALSE}
  ref<-readRDS(ref_path)
  if(missing(split_by)){
    data.seurat<-GEX
    query.projected <- make.projection(data.seurat, ref=ref, filter.cells = filtering)
    #Set up a list and not a dataframe if the whole GEX is used
    list_data<-query.projected
    query.projected<-list()
    query.projected[[1]]<-list_data
    #Projection
    split_by_character_vector<-c("Wohle VGM")
    projection_plots<-list()
    for(i in 1:length(query.projected)){
      projection_plots[[i]]<-print(plot.projection(ref,query.projected[[i]]) + ggtitle(split_by_character_vector[[i]]))
    }
    #Create two new columns in GEX data frame
    GEX<-VGM$GEX
    meta_data_new<-list()
    for (i in 1:length(query.projected)) {
      tils_prediction <- cellstate.predict(ref=ref, query=query.projected[[i]], )
      meta_data<-merge(GEX@meta.data, select(tils_prediction@meta.data, orig_barcode, functional.cluster,functional.cluster.conf), by = "orig_barcode")
      meta_data_new<-rbind(meta_data_new, meta_data)
    }
    #Recreate total GEX with functional.cluster
    GEX_total<-GEX@meta.data$orig_barcode
    GEX_barcode<-meta_data_new$orig_barcode
    out_filtered_barcode<-as.data.frame(setdiff(GEX_total, GEX_barcode))
    colnames(out_filtered_barcode)<-"orig_barcode"
    #extract outfilter row from vgm$GEX
    out_filtered_dataframe<-merge(GEX@meta.data, out_filtered_barcode, by = "orig_barcode")
    out_filtered_dataframe$functional.cluster <- NA
    out_filtered_dataframe$functional.cluster.conf <- NA
    GEX@meta.data<-rbind(meta_data_new, out_filtered_dataframe)
    #Barplots
    #With NA cells
    if(NA_cells == TRUE){
      cells<-GEX@meta.data
      cells$functional.cluster <- replace(cells$functional.cluster,is.na(cells$functional.cluster),"Undetermined")
      
    } else if (NA_cells == FALSE){
      cells<-GEX@meta.data[!is.na(GEX@meta.data$functional.cluster) | !is.na(GEX@meta.data$functional.cluster.conf), ]
    }
    barplots<-ggplot(cells, aes(x=split_by_character_vector, fill = functional.cluster)) + geom_bar(stat = "count", position = "fill", color="black", width = 0.6)  + theme_classic() +
      ylab("Fraction of cells") + xlab("Group") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ggtitle("") +
      scale_fill_manual(values = c("CD8_EarlyActiv" = "#F8766D", "CD8_EffectorMemory" = "#53B400", "CD8_NaiveLike" = "#00B6EB","CD8_Tex"="#edbe2a", "CD8_Tpex" = "#A58AFF", "Tfh" ="#FF0000", "Th1"="#87f6a5", "Treg"="#e812dd", "CD4_NaiveLike" = "#d1cfcc", "Undetermined" = "#000000")) + 
      scale_y_continuous(expand = c(0,0)) + labs(fill="predicted cell state") + scale_x_discrete(labels = split_by_character_vector)
    outputs<-list()
    outputs[[1]]<-GEX
    outputs[[2]]<-query.projected
    outputs[[3]]<-projection_plots
    outputs[[4]]<-barplots
    stop(return(outputs))
  }
  #Split by parameter
  subgroups<-unique(GEX@meta.data[[split_by]])
  subgroups<-na.omit(subgroups)
  if(length(subgroups)<=1){
    data.seurat<-GEX
    query.projected <- make.projection(data.seurat, ref=ref, filter.cells = filtering)
    #Set up a list and not a dataframe if the whole GEX is used
    list_data<-query.projected
    query.projected<-list()
    query.projected[[1]]<-list_data
    #Projection
    split_by_character_vector<-c("Wohle VGM")
    projection_plots<-list()
    for(i in 1:length(query.projected)){
      projection_plots[[i]]<-print(plot.projection(ref,query.projected[[i]]) + ggtitle(split_by_character_vector[[i]]))
    }
    #Create two new columns in GEX data frame
    GEX<-VGM$GEX
    meta_data_new<-list()
    for (i in 1:length(query.projected)) {
      tils_prediction <- cellstate.predict(ref=ref, query=query.projected[[i]], )
      meta_data<-merge(GEX@meta.data, select(tils_prediction@meta.data, orig_barcode, functional.cluster,functional.cluster.conf), by = "orig_barcode")
      meta_data_new<-rbind(meta_data_new, meta_data)
    }
    #Recreate total GEX with functional.cluster
    GEX_total<-GEX@meta.data$orig_barcode
    GEX_barcode<-meta_data_new$orig_barcode
    out_filtered_barcode<-as.data.frame(setdiff(GEX_total, GEX_barcode))
    colnames(out_filtered_barcode)<-"orig_barcode"
    #extract outfilter row from vgm$GEX
    out_filtered_dataframe<-merge(GEX@meta.data, out_filtered_barcode, by = "orig_barcode")
    out_filtered_dataframe$functional.cluster <- NA
    out_filtered_dataframe$functional.cluster.conf <- NA
    GEX@meta.data<-rbind(meta_data_new, out_filtered_dataframe)
    #Barplots
    #With NA cells
    if(NA_cells == TRUE){
      cells<-GEX@meta.data
      cells$functional.cluster <- replace(cells$functional.cluster,is.na(cells$functional.cluster),"Undetermined")
      
    } else if (NA_cells == FALSE){
      cells<-GEX@meta.data[!is.na(GEX@meta.data$functional.cluster) | !is.na(GEX@meta.data$functional.cluster.conf), ]
    }
    barplots<-ggplot(cells, aes(x=split_by_character_vector, fill = functional.cluster)) + geom_bar(stat = "count", position = "fill", color="black", width = 0.6)  + theme_classic() +
      ylab("Fraction of cells") + xlab("Group") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ggtitle("") +
      scale_fill_manual(values = c("CD8_EarlyActiv" = "#F8766D", "CD8_EffectorMemory" = "#53B400", "CD8_NaiveLike" = "#00B6EB","CD8_Tex"="#edbe2a", "CD8_Tpex" = "#A58AFF", "Tfh" ="#FF0000", "Th1"="#87f6a5", "Treg"="#e812dd", "CD4_NaiveLike" = "#d1cfcc", "Undetermined" = "#000000")) + 
      scale_y_continuous(expand = c(0,0)) + labs(fill="predicted cell state") + scale_x_discrete(labels = split_by_character_vector)
    outputs<-list()
    outputs[[1]]<-GEX
    outputs[[2]]<-query.projected
    outputs[[3]]<-projection_plots
    outputs[[4]]<-barplots
    stop(return(outputs))
  }
  data.seurat.list<-SplitObject(GEX, split.by = split_by)
  #Run projecTILS algorithm, this step takes really time
  query.projected.list <- make.projection(data.seurat.list, ref=ref, filter.cells = filtering)
  #Set up a list and not a dataframe if the whole GEX is used
  if(length(query.projected.list)==1){
    list_data<-query.projected.list
    query.projected.list<-list()
    query.projected.list[[1]]<-list_data
  }
  #Projection
  split_by_character_vector<-names(query.projected.list)
  projection_plots<-list()
  for(i in 1:length(query.projected.list)){
    projection_plots[[i]]<-print(plot.projection(ref,query.projected.list[[i]]) + ggtitle(split_by_character_vector[[i]]))
  }
  #Barplots fraction of cells with predicted state per cluster
  #Create two new columns in GEX data frame
  GEX<-VGM$GEX
  meta_data_new<-list()
  for (i in 1:length(query.projected.list)) {
    tils_prediction <- cellstate.predict(ref=ref, query=query.projected.list[[i]], )
    meta_data<-merge(GEX@meta.data, select(tils_prediction@meta.data, orig_barcode, functional.cluster,functional.cluster.conf), by = "orig_barcode")
    meta_data_new<-rbind(meta_data_new, meta_data)
  }
  #Recreate total GEX with functional.cluster
  GEX_total<-GEX@meta.data$orig_barcode
  GEX_barcode<-meta_data_new$orig_barcode
  out_filtered_barcode<-as.data.frame(setdiff(GEX_total, GEX_barcode))
  colnames(out_filtered_barcode)<-"orig_barcode"
  #extract outfilter row from vgm$GEX
  out_filtered_dataframe<-merge(GEX@meta.data, out_filtered_barcode, by = "orig_barcode")
  out_filtered_dataframe$functional.cluster <- NA
  out_filtered_dataframe$functional.cluster.conf <- NA
  GEX@meta.data<-rbind(meta_data_new, out_filtered_dataframe)
  #Barplots
  #With NA cells
  if(NA_cells == TRUE){
    cells<-GEX@meta.data
    cells$functional.cluster <- replace(cells$functional.cluster,is.na(cells$functional.cluster),"Undetermined")
    
  } else if (NA_cells == FALSE){
    cells<-GEX@meta.data[!is.na(GEX@meta.data$functional.cluster) | !is.na(GEX@meta.data$functional.cluster.conf), ]
  }
  barplots<-ggplot(cells, aes(x = as.character(.data[[split_by]]) , fill = functional.cluster))+ geom_bar(stat = "count", position = "fill", color="black", width = 0.6)  + theme_classic() +
    ylab("Fraction of cells") + xlab(split_by) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ggtitle("") +
    scale_fill_manual(values = c("CD8_EarlyActiv" = "#F8766D", "CD8_EffectorMemory" = "#53B400", "CD8_NaiveLike" = "#00B6EB", 
                                 "CD8_Tex"="#edbe2a", "CD8_Tpex" = "#A58AFF", "Tfh" ="#FF0000", "Th1"="#87f6a5", "Treg"="#e812dd", "CD4_NaiveLike" = "#d1cfcc", "Undetermined" = "#000000")) + 
    scale_y_continuous(expand = c(0,0)) + labs(fill="predicted cell state") + scale_x_discrete(labels = c(split_by_character_vector))
  #Outputs
  outputs<-list()
  outputs[[1]]<-GEX
  outputs[[2]]<-query.projected.list
  outputs[[3]]<-projection_plots
  outputs[[4]]<-barplots
  return(outputs)
}
