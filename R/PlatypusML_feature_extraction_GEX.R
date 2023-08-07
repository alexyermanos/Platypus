

#' This PlatypusML_feature_extraction_GEX function takes as input specified features from the second output of the VDJ_GEX_matrix function and encodes
#' according to the specified strategy. The function returns a matrix containing the encoded extracted features as columns and the different cells as rows.
#' This function should be called as a first step in the process of modeling the VGM data using machine learning.
#' @param VGM output of the VDJ_GEX_matrix function, containing both VDJ and GEX objects.
#' @param encoding.level String. Specifies on which level the features will be extracted. There are three possible options: "clone" (one random sample per clone),
#' "clone.avg" (average expression per clone), "unique.sequence" (selecting only unique sequences based on a specified sequence (in the unique.sequence argument)).
#' Defaults to "clone.avg".
#' @param unique.sequence String. Needs to be specified only when encoding.level is set to "unique.sequence". The name of the sequence on which unique selection should be based on.
#' Defaults to "VDJ_cdr3s_aa".
#' @param which.features String. Information on which GEX features should be encoded. Options are "varFeatures" (the 1000 most variable features obtained by Seurat::FindVariableFeatures)
#' or "PCs" (the top n PCs, number of PCs to be defined in n.PCs). Defaults to "PCs".
#' @param n.PCs Integer. Number of PCs to be used if choosing which.features == "PCs". Max 50. Defaults to 20.
#' @param which.label String. The name of the column in VGM[[2]] which will be appended to the encodings and used as a label in a chosen ML model later.
#' The label has to be a binary label. If missing, no label will be appended to the encoded features.
#' @param problem String ("classification" or "regression"). Whether the return matrix will be used in a classification problem or a regression one. Defaults to "classification".
#' @param verbose.classes Boolean. Whether to display information on the distribution of samples between classes. Defaults to TRUE.
#' For this parameter to be set to TRUE, classification must all be set to TRUE (default).
#' @param platypus.version This function works with "v3" only, there is no need to set this parameter.
#' @return A dataframe containing the encoded features and its label, each row corresponding to a different cell.
#' The label can be found in the last column of the dataframe returned. If which.label="NA" only the encoded features are returned.
#' @export
#' @examples
#' \dontrun{
#' To return the encoded gene expression in form of
#' the 20 PCs at the clone level
#' (average expression per clone).
#' Attaching the "GP33_binder" label to be used in downstream ML models.
#' features_PCs_GP33_binder <- PlatypusML_feature_extraction_GEX(
#' VGM = VGM,
#' encoding.level = "clone.avg",
#' which.features = "PCs",
#' n.PCs = 20,
#' which.label = "GP33_binder")
#'}


PlatypusML_feature_extraction_GEX <- function(VGM,
                                          encoding.level,
                                          unique.sequence,
                                          which.features,
                                          n.PCs,
                                          which.label,
                                          problem,
                                          verbose.classes,
                                          platypus.version){


  ############### checking for missing/invalid inputs & setting defaults ##############################
  unique_sequence <- NULL
  clonotype_distinct <- NULL
  bind_rows <- NULL
  if (missing(VGM)) stop("Please provide VGM input for this function")

  #encoding.level
  if (missing(encoding.level)){
    encoding.level <- "clone.avg"
    print("encoding.level was set to clone.avg")}
  if (encoding.level == "unique.sequence" & missing(unique.sequence)){
    unique.sequence <- "VDJ_cdr3s_aa"
    print("unique.sequence was set to VDJ_cdr3s_aa")}
  if (missing(unique.sequence)) unique.sequence <- "NA"
  if ( unique.sequence != "NA" & !(unique.sequence %in% names(VGM[[2]]@meta.data))) stop("Please provide a valid column (unique.sequence) that exists in VGM[[2]]@meta.data")
  if (!(encoding.level %in% c("clone", "clone.avg", "unique.sequence", "cell"))){
    stop("chosen encoding.level not found. Please provide one of the following options: clone, clone.avg, unique.sequence", "cell")}

  #classification & which.label
  if (missing(problem)) problem <- "classification"
  if (problem == "classification" & missing(which.label)) stop("Plase provide a label (which.label) for this function")
  if (problem !="classification"  & missing(which.label)) which.label <- "NA"
  if (!(which.label %in% names(VGM[[2]]@meta.data))) stop("Please provide a valid label (which.label) that exists in VGM[[2]]@meta.data")
  if (missing(verbose.classes)) verbose.classes <- TRUE

  #which.features
  if (!(which.features %in% c("PCs", "varFeatures"))){
    stop ("chosen which.features not found. Please provide one of the following options: PCs, varFeatures")
  }
  if (missing(which.features)){
    which.features <- "PCs"
    print("which.features was set to PCs")}
  if (missing(n.PCs) & which.features=="PCs"){
    n.PCs <- 20
    print("n.PCs was set to 20")}

  ############### encoding.level == "clone" OR "unique.sequence" OR "cell" ##############################
  if (encoding.level == "unique.sequence" | encoding.level == "clone" | encoding.level == "cell"){

    #filter the cells based on the selected encoding method

    if (encoding.level == "unique.sequence"){
      VGM[[2]]$unique_sequence <- VGM[[2]]@meta.data[unique.sequence]
      unique_sequences <- unique(VGM[[2]]$unique_sequence)
      unique_sequences_idx <- as.vector(match(unique_sequences,VGM[[2]]$unique_sequence))
      VGM[[2]]$unique_sequence <-VGM[[2]]$unique_sequence[unique_sequences_idx]
      #subset Seurat object (GEX)
      VGM[[2]] <- subset(VGM[[2]], subset=unique_sequence != "NA")
    }

    if (encoding.level == "clone"){
      VGM[[2]]$clonotype_distinct <- paste(VGM[[2]]$clonotype_id_10x, VGM[[2]]$sample_id)
      unique_clones <- unique(VGM[[2]]$clonotype_distinct)
      #remove rows with clonotype_id_10x==NA
      unique_clones <- unique_clones[!(grepl("NA", unique_clones))]
      unique_clones_idx <- as.vector(match(unique_clones,VGM[[2]]$clonotype_distinct))
      VGM[[2]]$clonotype_distinct <- VGM[[2]]$clonotype_distinct[unique_clones_idx]
      #subset Seurat object (GEX)
      VGM[[2]] <- subset(VGM[[2]], subset=clonotype_distinct != "NA")
    }

    #if the encoding level is "cell" no selection is required and the encoding will be done for all of the sequences

    #remove cells with label==NA
    VGM[[2]]$label <- VGM[[2]]@meta.data[which.label]
    #subset Seurat Object (GEX)
    VGM[[2]] <- subset(VGM[[2]], subset=label != "NA")


    #class balance summary
    if (problem == "classification" & verbose.classes){
      labels <- VGM[[2]]@meta.data[which.label]
      classes <-  unique(labels)
      for (i in 1:nrow(classes)){
        count_class <- length(which(labels == classes[i,]))
        perc_class <- round(count_class/nrow(labels)*100,2)
        print(sprintf("Class %s contains %s samples, representing %s%% of the total samples.", classes[i,], count_class, perc_class))
      }
    }

    # PCs encoding (selecting the top n.PCs PCs)
    if (which.features=="PCs"){
      VGM[[2]] <- Seurat::RunPCA(VGM[[2]], features=SeuratObject::VariableFeatures(object = VGM[[2]]), verbose=FALSE)
      encoding <- VGM[["GEX"]]@reductions[["pca"]]@cell.embeddings[,1:n.PCs]
    }

    # varFeatures encoding (selecting the top 1000 most variable genes)
    if (which.features=="varFeatures"){
      VGM[[2]] <- Seurat::FindVariableFeatures(VGM[[2]], nfeatures=1000, verbose=FALSE)
      rownames <- rownames(VGM[[2]]@assays[["RNA"]]@counts)
      varfeatures <- VGM[[2]]@assays[["RNA"]]@var.features
      varfeatures <- rownames %in% varfeatures
      encoding <- VGM[[2]]@assays[["RNA"]]@counts[varfeatures,]
      encoding <- as.data.frame(t(encoding))
    }

    if (which.label!="NA"){
      # append label to encoding
      label <- VGM[[2]]@meta.data[which.label]
      encoding <- cbind(encoding, label)
      encoding <- as.data.frame(encoding)
      return (encoding)
    } else {
      # return encoding only (without an appended label)
      return (as.data.frame(encoding))
    }


    ############### encoding.level == "clone.avg" ##############################
  } else {

    #filter the cells based on the selected encoding method "clone.avg"

    #remove cells with label==NA
    VGM[[2]]@meta.data$label <- VGM[[2]]@meta.data[which.label]
    #subset Seurat object (GEX)
    VGM[[2]] <- subset(VGM[[2]], subset=label != "NA")

    VGM[[2]]$clonotype_distinct <- paste(VGM[[2]]$clonotype_id_10x, VGM[[2]]$sample_id)
    #remove rows with clonotype_id_10x==NA
    VGM[[2]]$clonotype_distinct[grepl("NA",VGM[[2]]$clonotype_distinct)] <- NA
    #subset Seurat object (GEX)
    VGM[[2]] <- subset(VGM[[2]], subset=clonotype_distinct != "NA")
    #find out which samples had which binary label
    label_0 <- as.vector(unique(VGM[[2]]$sample_id[VGM[[2]]@meta.data$label==0]))
    label_1 <- as.vector(unique(VGM[[2]]$sample_id[VGM[[2]]@meta.data$label==1]))
    sample_ids <- as.vector(unique(VGM[[2]]$sample_id))
    #average expression per clone, subsetting SeuratObject (GEX)
    VGM[[2]] <- Seurat::AverageExpression(object=VGM[[2]],
                                          return.seurat=TRUE,
                                          group.by = "clonotype_distinct",
    )
    #add labels to new Seurat object
    VGM[[2]]@meta.data$label <- NA
    for (i in label_0){
      VGM[[2]]@meta.data$label[grepl(i,names(VGM[[2]]@active.ident))] <- 0
    }
    for (j in label_1){
      VGM[[2]]@meta.data$label[grepl(j,names(VGM[[2]]@active.ident))] <- 1
    }

    #class balance summary
    if (problem == "classification" & verbose.classes){
      labels <- VGM[[2]]@meta.data$label
      classes <-  unique(labels)
      for (i in 1:length(classes)){
        count_class <- length(which(labels == classes[i]))
        perc_class <- round(count_class/length(labels)*100,2)
        print(sprintf("Class %s contains %s samples, representing %s%% of the total samples.", classes[i], count_class, perc_class ))
      }

    }

    # PCs encoding (selecting the top n.PCs PCs)
    if (which.features=="PCs"){
      VGM[[2]] <- Seurat::FindVariableFeatures(VGM[[2]], verbose=FALSE)
      VGM[[2]] <- Seurat::RunPCA(VGM[[2]], features=SeuratObject::VariableFeatures(object = VGM[[2]]), verbose=FALSE)
      encoding <- VGM[["GEX"]]@reductions[["pca"]]@cell.embeddings[,1:n.PCs]
    }

    # varFeatures enoding (selecting the top 1000 most variable genes)
    if (which.features=="varFeatures"){
      VGM[[2]] <- Seurat::FindVariableFeatures(VGM[[2]], nfeatures=1000, verbose=FALSE)
      rownames <- rownames(VGM[[2]]@assays[["RNA"]]@counts)
      varfeatures <- VGM[[2]]@assays[["RNA"]]@var.features
      varfeatures <- rownames %in% varfeatures
      encoding <- VGM[[2]]@assays[["RNA"]]@counts[varfeatures,]
      encoding <- as.data.frame(t(encoding))
    }

    if (which.label!="NA"){
      # append label to encoding
      label_col <- VGM[[2]]@meta.data$label
      colnames_encoding <- colnames(encoding)
      encoding <- cbind(encoding, label_col)
      encoding <- as.data.frame(encoding)
      colnames(encoding) <- c(colnames_encoding, which.label)
      return (encoding)
    } else {
      # return encoding only (without an appended label)
      return (as.data.frame(encoding))
    }
  }

}
