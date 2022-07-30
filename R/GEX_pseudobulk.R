#' Function that performs pseudo-bulking on the data (VGM input), according to criteria specified by the User, and uses the pseudo-bulked data to perform Differential Gene Expression (DGE) analysis.
#' @param vgm.input Output of the VDJ_GEX_matrix function. Mandatory
#' @param column.group Character vector. Mandatory. Column name of VDJ_GEX_matrix[[2]] where the groups to be tested for differenetial gene expression are located
#' @param group1 Strings vector. Mandatory. Samples to be grouped together for differential expression analysis against group2 (if pool=TRUE or column.comparison!=NULL). If pool=FALSE, vector containing samples to be tested individually against samples with the same index in the vector of group2.
#' @param group2 Strings vector. Mandatory. Samples to be grouped together for differential expression analysis against group1 (if pool=TRUE or column.comparison!=NULL). If pool=FALSE, vector containing samples to be tested individually against samples with the same index in the vector of group1.
#' @param column.comparison Character vector. Defaults to NULL. Column name of VDJ_GEX_matrix[[2]] where the comparison cathegories are located, if DGE between group1 and group2 is performed across diferrent cathegories.
#' @param comparison, Strings vector. Defaults to NULL. Comparison cathegories, if more than one is present.
#' @param pool Logical. Defaults to FALSE. Indicates whether samples specified in group1 and group2 are to be pooled together within the same group.
#' @param platypus.version This function works with "v3" only, there is no need to set this parameter
#' @return A data.frame or list of data.frames containing the results of the DGE analysis for every level of pseudo-bulking.
#' @export
#' @examples
#' \dontrun{
#'pseudo_DE<-GEX_pseudobulk(
#'vgm.input=VGM_NP396_GP33_GP66_labelled,
#'column.group="sample_id",
#'group1 = c("s1"), group2 = c("s4"),
#'column.comparison = "seurat_clusters",
#'comparison=c(2,3,5,9), pool=FALSE)
#'}
#'

GEX_pseudobulk<-function(vgm.input, #mandatory  #the function takes vgm as an input
                         column.group, #column name of vgm[[2]] where the groups to be tested for differenetial gene expression are located
                         group1, #mandatory #samples to be grouped together for differential expression analysis against group2
                         group2, #mandatory #samples to be grouped together for differential expression analysis against group1
                         column.comparison, #column name of vgm[[2]] where the comparison cathegories are located, if DGE is performed across diferrent cathegories (eg sample_id as column.group and seurat_clusters as column.comparison). By default set to NULL
                         comparison, #vector indicating the comparison cathegories, if more than one is present. By default set to NULL
                         pool, #logical indicating whether the samples in group1 and group2 are to be pooled together within the same group. By default FALSE
                         platypus.version
){

  platypus.version <- "v3"

  #For RMD check
  . <- NULL
  PValue <- NULL
  FDR <- NULL
  p_adj <- NULL
  logFC <- NULL

  #set default value for column.group
  if(missing(column.group)){
    column.group="seurat_clusters"
  }

  #set default value for column.comparison
  if(missing(column.comparison)){
    column.comparison=NULL
  }

  #set default value for comparison
  if(missing(comparison)){
    comparison=NULL
  }

  #set default value for pool
  if(missing(pool)){
    pool=FALSE
  }

  groups_together<-c(group1, group2)

  #take subset of vgm[[2]] according to the input column.group and input groups given by the User

  where_groups<-c()
  for(i in groups_together){
    #get the index of the groups
    where_groups<-append(where_groups, which(vgm.input[[2]][[column.group]]==as.character(i)))
  }

  input<-vgm.input[[2]][,sort(where_groups)]

  #if theere are more samples in each comparison group, create a new column in the meta.data and pull them together according to the grouping (under the name "group1" and "group2" respectively)

  input@meta.data$groups<-rep(0,length(input[[column.group]]))

  for(i in group1){
    temp<-which(input[[column.group]]==as.character(i))
    input@meta.data$groups[temp]<-"group1"
  }

  for(i in group2){
    temp<-which(input[[column.group]]==as.character(i))
    input@meta.data$groups[temp]<-"group2"
  }

  input$groups<-factor(input$groups)


  #if column.comparison and comparison are NOT NULL, select only data corresponding to the given cathegories and indicate them in the new column. These will represent the levels of comparison between the two groups . Useful when one wants to perform DGE analysis according to different cathegories. column.group and column.comparison need to be different to avoid overlappings.
  #in case there are multiple samples per group, these are pooled together by being asssigned to group1 or group2

  if(!(is.null(column.comparison))&!(is.null(comparison))& pool==FALSE){

    #create column containing comparison levels
    input@meta.data$comparison<-rep(0,length(input[[column.comparison]]))
    where_comparison<-c()
    for(i in comparison){
      #get the index of the groups
      temp<-c()
      temp<-which(input[[column.comparison]]==as.character(i))
      where_comparison<-append(where_comparison, temp)
      input@meta.data$comparison[temp]<-paste(c("comparison", i), collapse = "_")
    }
    input$comparison<- factor(input$comparison)
    input<-input[,sort(where_comparison)]

    #create a column containing the comparison-group combination
    levels=c()
    for(g in 1:length(input$groups)){
      levels<-append(levels,paste( c(input$groups[g], input$comparison[g]), collapse="_"))
    }
    input$levels<-levels

    #input object for summarizeAssayByGroup() needs to be a Summarised Experiment
    input<-Seurat::as.SingleCellExperiment(input)
    assign("input", input, globalenv())

    # aggregate by cluster-sample
    g <- input@colData[, c("comparison", "groups")]
    pb <- Matrix.utils::aggregate.Matrix(t(input@assays@data@listData$counts),
                           groupings = g, fun = "sum")

    #build the vector defining the grouping for split.data.frame() function
    f<-c()
    for(i in unique(input$comparison)){
      f<-append(f, rep(i, length(levels(input$groups))))
    }
    # split by cluster, transform & rename columns
    pb <- split.data.frame(pb, f)

    for(i in unique(input$comparison)){
      pb[[i]] <- magrittr::set_colnames(t(pb[[i]]), unname(levels(input$groups)))
    }
    #construct a SCE where each assay sheet corresponds to one comparison level: assays = comparison, rows = genes, columns = groups
    pb <- SingleCellExperiment::SingleCellExperiment(assays = pb)

    #correspondences between comparison, levels and groups id
    df_summary<-unique(input@colData[,c(as.character(column.group), "comparison", "groups", "levels")])

    df_analysis<-unique(input@colData[,c( "comparison", "groups", "levels")])

    #construct design matrix
    design <- stats::model.matrix(~ 0 + df_analysis$groups)
    design <-magrittr::set_rownames(design, df_analysis$levels)
    design <-magrittr::set_colnames(design, levels(df_analysis$groups))
    contrast <- limma::makeContrasts("group2-group1", levels = design)


    res <- lapply(as.character(unique(input$comparison)), function(k) {

      . <- NULL
      PValue <- NULL
      FDR <- NULL
      p_adj <- NULL
      logFC <- NULL

      y<-pb@assays@data@listData[[k]]
      y<-edgeR::DGEList(y, remove.zeros = TRUE)
      y <- edgeR::calcNormFactors(y)

      #dispersion cannot be calculted without replicate, so if only two pseudo-bulked groups are compared, dispersion needs to be estimated.Ideally you would have biologically replicates though.
      y <- edgeR::estimateGLMCommonDisp(y, method="deviance", robust=TRUE, subset=NULL)
      #subset design matrix according to comparison level
      fit <- edgeR::glmFit(y, design[which(rownames(design)==df_analysis$levels[df_analysis$comparison==k][1]|rownames(design)==df_analysis$levels[df_analysis$comparison==k][2]),])
      fit <- edgeR::glmLRT(fit, contrast = contrast)
      edgeR::topTags(fit, n = Inf, sort.by = "none")$table %>%
        dplyr::mutate(gene = rownames(.), comparison_id = k) %>%
        dplyr::rename(p_val = PValue, p_adj = FDR) %>% dplyr::filter(p_adj < 0.1, abs(logFC) > 1) %>% dplyr::arrange(p_adj)
    })
    # number and % of differential genes by comparison cluster
    number_de <- vapply(res, nrow, numeric(1))
    res$number_de<-cbind(number_de, percentage_de = number_de / nrow(input) * 100)
    #create a df with correspondences between column.group, comparison and groups id, for easier visualization
    df_summary<-unique(input@colData[,c(as.character(column.group), "comparison", "groups", "levels")])
    res$summary<-df_summary

    return(res)

    #if groups contain each more than one sample and pool is set to TRUE, these are pooled together according to their belonging to group1 or group2.
  }else if((length(group1)>1|length(group2)>1) & pool==TRUE){

    #input object for summarizeAssayByGroup() needs to be a Summarised Experiment
    input<-Seurat::as.SingleCellExperiment(input)

    #sum counts across all cells in each cluster to obtain “pseudo-bulk” samples for differential expression analyses between conditions.
    pseudo<-scuttle::aggregateAcrossCells(input, ids=input@colData[,c(as.character(column.group))])

    #construct design matrix
    design <- stats::model.matrix(~ 0 + pseudo$groups)
    design<-magrittr::set_rownames(design, pseudo[[column.group]])
    design<-magrittr::set_colnames(design, levels(pseudo$groups))
    contrast <- limma::makeContrasts("group2-group1", levels = design)

    #perform DGE analysis
    y<-pseudo@assays@data@listData$counts
    y<-edgeR::DGEList(y, remove.zeros = TRUE)
    y <- edgeR::calcNormFactors(y)
    #dispersion cannot be calculted without replicate, so if only two pseudo-bulked groups are compared, dispersion needs to be estimated
    y <- edgeR::estimateGLMCommonDisp(y, method="deviance", robust=TRUE, subset=NULL)
    #subset design matrix according to comparison level
    fit <- edgeR::glmFit(y, design)
    fit <- edgeR::glmLRT(fit, contrast = contrast)
    toptags<-edgeR::topTags(fit, n = Inf, sort.by = "none")$table %>%
      dplyr::mutate(gene = rownames(.)) %>%
      dplyr::rename(p_val = PValue, p_adj = FDR) %>% dplyr::filter(p_adj < 0.05, abs(logFC) > 1) %>% dplyr::arrange(p_adj)

    return(toptags)

    #if pool is set to FALSE and the amount of samples in group1 and group2 is equal, perform DGE analysis comparing the first sample in group1 with the first one in group2, the second sample in group1 with the second one in group2 and so on... group1 amd group2 need to contain different and non repeated samples
  }else if(pool==FALSE & length(group1)==length(group2)){

    #create a column containing the comparison-group combination
    input$comparison=rep(0, length(input[[column.group]]))
    for(g in 1:length(group1)){
      temp1<-which(input[[column.group]]==as.character(group1[g]))
      temp2<-which(input[[column.group]]==as.character(group2[g]))
      input$comparison[temp1]<-paste(c("comparison", g), collapse = "_")
      input$comparison[temp2]<-paste(c("comparison", g), collapse = "_")
    }

    input$comparison<-factor(input$comparison)

    #input object for summarizeAssayByGroup() needs to be a Summarised Experiment
    input<-Seurat::as.SingleCellExperiment(input)

    # aggregate by cluster-sample
    g <- input@colData[, c("comparison", as.character(column.group))]
    pb <-  Matrix.utils::aggregate.Matrix(t(input@assays@data@listData$counts),
                           groupings = g, fun = "sum")

    # split by cluster, transform & rename columns
    pb <- split.data.frame(pb, rep(unique(input$comparison), length(levels(input$groups)))) %>%
      lapply(function(u) magrittr::set_colnames(t(u), unname(levels(input$groups))))
    #construct a SCE where each assay sheet corresponds to one comparison level: assays = comparison, rows = genes, columns = groups
    pb <- SingleCellExperiment::SingleCellExperiment(assays = pb)

    #create data.frame with the correspondences between column.group, comparison and groups id, for easier visualization
    df_summary<-unique(input@colData[,c(as.character(column.group), "comparison", "groups")])

    #construct design matrix
    design <- stats::model.matrix(~ 0 + df_summary$groups)
    design<-magrittr::set_rownames(design, df_summary[[column.group]])
    design<-magrittr::set_colnames(design, levels(df_summary$groups))
    contrast <- limma::makeContrasts("group2-group1", levels = design)

    #perform DGE analysis
    res <- lapply(levels(input$comparison), function(k) {

      . <- NULL
      PValue <- NULL
      FDR <- NULL
      p_adj <- NULL
      logFC <- NULL

      y<-pb@assays@data@listData[[k]]
      y<-edgeR::DGEList(y, remove.zeros = TRUE)
      y <- edgeR::calcNormFactors(y)
      #dispersion cannot be calculted without replicate, so if only two pseudo-bulked groups are compared, dispersion needs to be estimated
      y <- edgeR::estimateGLMCommonDisp(y, method="deviance", robust=TRUE, subset=NULL)
      #subset design matrix according to comparison level
      fit <- edgeR::glmFit(y, design[which(rownames(design)==df_summary[[column.group]][df_summary$comparison==k][1]|rownames(design)==df_summary[[column.group]][df_summary$comparison==k][2]),])
      fit <- edgeR::glmLRT(fit, contrast = contrast)
      edgeR::topTags(fit, n = Inf, sort.by = "none")$table %>%
        dplyr::mutate(gene = rownames(.), comparison_id = k) %>%
        dplyr::rename(p_val = PValue, p_adj = FDR) %>% dplyr::filter(p_adj < 0.05, abs(logFC) > 1) %>% dplyr::arrange(p_adj)
    })
    # number and % of differential genes by comparison cluster
    number_de <- vapply(res, nrow, numeric(1))
    res$number_de<-cbind(number_de, percentage_de = number_de / nrow(input) * 100)
    res$summary<-df_summary

    return(res)

  }else{
    return("Error: at least two samples per group are needed, if no comparison levels are provided. If pull=FALSE, the lengths of group1 and group2 need to be equal.")
  }

}
