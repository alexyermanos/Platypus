#' Utility for feature barcode assignment including clonal information
#'
#'@description The VGM_expand_featurebarcodes function function can be used to trace back the cell origin of each sample after using cell hashing for single-cell sequencing. Replaces the original sample_id column of a vgm object with a pasted version of the original sample_id and the last digits of the feature barcode.
#'
#'The original sample_id is stored in a new column called original_sample_id. Additionally, a second new column is created containing final barcode assignment information.
#'Those barcodes match the origin FB_assignment if by.majority.barcodes is set to FALSE (default). However, if this input parameter is set to TRUE, the majority barcode assignment in stored in this colum.
#'
#'Note: The majority barcode of a cell is the feature barcode which is  most frequently assigned to the cells clonotype (10x default clonotype).
#'The majority barcode assignment can be used under the assumption that all cells which are assigned to the same clonotype (within one sample), originate from the same donor organ or at least the same donor depending on the experimental setup.
#'
#'For example: The original sample_id of a cell is "s1", the cell belongs to "clonotype1" and the feature barcode assigned to it is "i1-TotalSeq-C0953". If by.majority.barcodes default (FALSE) is used, the resulting new sample_id would be "s1_0953".
#'However, if majority barcode assignment is used AND "i1-TotalSeq-C0953" is not the most frequently occurring barcode in "clonotype1" but rather barcode "i1-TotalSeq-C0951", the new sample_id would be "s1_0951".
#'--> e.g., if 15 cells belong to clonotype1: 3 cells have no assigned barcode, 2 are assigned to "i1-TotalSeq-C0953" and 10 are assigned to "i1-TotalSeq-C0951" --> all 15 cells will have the new sample_id "s1_0951".
#'@param vgm VGM output of VDJ_GEX_matrix function (Platypus V3)
#'@param by.majority.barcodes Logical. Default is FALSE. Indicated whether strict barcode assignment or majority barcode assignment should be used to create the new sample_id. If TRUE, for each clonotype the most frequent feature barcode will be chosen and assigned to each cell, even if that cell itself does not have this particular barcode assigned.
#'@param integrate.in.gex Logical. Default is FALSE. If TRUE, the newly created sample_id's are integrated into gex component as well. Not recommended if no further gex analysis is done due to much longer computational time.
#'@param vdj.only Logical. Defines if only vdj information is provided as input. Default is set to FALSE. If set to TRUE a vdj dataframe has to be provided as input (vgm = vdj). Also, integrate.in.gex is automatically set to FALSE since no gex (vgm[[2]]) information is provided.
#'@param platypus.version This function works with "v3" only, there is no need to set this parameter.
#'@return This function returns a vgm with new sample_id's in case vdj.only is set to FALSE (default). If vdj.only is set to true only the vdj dataframe with new sample_id's is returned.
#'Note: If vdj.only is set to default (FALSE), VDJ information in the metadata of the GEX object is necessary. For this set integrate.VDJ.to.GEX to TRUE in the VDJ_GEX_matrix function
#'
#'@export
#'@examples
#' #For Platypus version 3
#'
#' # 1. If only vdj data (vgm[[1]]) and
#' #strict feature barcode assignment is used:
#' vgm_expanded_fb <- VGM_expand_featurebarcodes(
#' vgm = small_vgm[[1]],
#' by.majority.barcodes = FALSE,
#' integrate.in.gex=FALSE, vdj.only= TRUE)
#'
#' # 2. If whole vgm and strict fb assignment is used
#' #(gex and vdj - necessary if gene expression analysis
#' # of sub-samples is desired):
#' vgm_expanded_fb <- VGM_expand_featurebarcodes(
#' vgm = small_vgm,
#' by.majority.barcodes = FALSE,
#' integrate.in.gex=TRUE, vdj.only= FALSE)
#'
#' # 3. If whole vgm and majority barcode assignement is used
#' #(gex and vdj) - necessary if gene expression analysis
#' #of sub-samples is desired):
#' vgm_expanded_fb <- VGM_expand_featurebarcodes(vgm = small_vgm,
#' by.majority.barcodes = TRUE,
#' integrate.in.gex=TRUE, vdj.only= FALSE)
#'
#' #Note: Majority barcode assignment is recommended
#' #if the assumption that all cells within one clonotype
#' #originate from the same sample sub-group is feasible.


VGM_expand_featurebarcodes <- function(vgm,
                                         by.majority.barcodes,
                                         integrate.in.gex,
                                         vdj.only,
                                         platypus.version){

  if(missing(vgm)) stop("Please provide vgm input for this function")

  if(missing(by.majority.barcodes)) by.majority.barcodes <- FALSE
  if(missing(integrate.in.gex)) integrate.in.gex <- FALSE
  if(missing(vdj.only)) vdj.only <- FALSE
  if(missing(platypus.version)) platypus.version <- "v3"


  #extract vdj
  if(vdj.only == FALSE){
    VDJ <- vgm[[1]]
  }else{
    VDJ <- vgm
    integrate.in.gex <- FALSE
  }

  sample.names <- unique(VDJ$sample_id)

  #create original_sample_id column in vdj dataframe to store the old sample_id's
  original_sample_id <- NULL
  VDJ$original_sample_id <- VDJ$sample_id
  VDJ$sample_id <- as.character(VDJ$sample_id) #turn sample_id into a charcter vector to be able to add new sample_id names later

  #create final_fb_assignment column
  final_fb_assignment <- NULL
  VDJ$final_fb_assignment <- NA

  # Step 1: Barcode Assignment
  if(by.majority.barcodes == FALSE){ #original barcode assigned

    VDJ$final_fb_assignment <- VDJ$FB_assignment

  }else{ #majority barcode assignment

    #iterate over every sample, create a list of unique feature barcodes and clonotypes for each sample
    for (i in 1:length(sample.names)) {
      curr_sample <- subset(VDJ, VDJ$sample_id == sample.names[i])

      curr_feature_barcodes <- unique(curr_sample$FB_assignment)
      curr_feature_barcodes <- curr_feature_barcodes[order(curr_feature_barcodes)]
      curr_clonotypes <- unique(curr_sample$clonotype_id_10x)

      #iterate over every clonotype, extract dataframe subsets of each clone, find the most frequently assigned feature barcode and assign this fb to the final_fb_assignment column for every cell in the df
      for(j in 1:length(curr_clonotypes)){
        curr_clone <- subset(curr_sample, curr_sample$clonotype_id_10x == curr_clonotypes[j])

        fb_freq_table <-as.data.frame(table(curr_clone$FB_assignment))
        fb_freq_table <- fb_freq_table[order(-fb_freq_table$Freq),] #sort table by fb frequency

        most_frequent_fb <- NULL
        if(nrow(fb_freq_table)>=2 & fb_freq_table$Var1[1] == "Not assignable"){ #In case there are two feature barcodes with the same frequency in this clonotype and the top list entry is "Not assignable", use the second entry
          if(fb_freq_table$Freq[1] == fb_freq_table$Freq[2]){
            most_frequent_fb <- as.character(fb_freq_table$Var1[2])
          }else{
            most_frequent_fb <- as.character(fb_freq_table$Var1[1])
          }
        }else{
          most_frequent_fb <- as.character(fb_freq_table$Var1[1]) #If there is only one fb with the higest frequency, simply use the first element of the ordered frequency table
        }

        #assign the most_frequent_fb to the final_fb_assignment column of all cells within the current clonotype and sample
        VDJ$final_fb_assignment[VDJ$sample_id == sample.names[i] & VDJ$clonotype_id_10x == curr_clonotypes[j]] <- most_frequent_fb
      }
    }
  }


  #Step 2: Overwrite old sample_id with pasted sample_id + final_fb_assignment information

  #Note: only last 4 digits of the feature barcodes are used
  for (i in 1:nrow(VDJ)){
    short_fb <- substring(VDJ$final_fb_assignment[i], nchar(VDJ$final_fb_assignment[i])-3)
    if(short_fb == "able"){
      short_fb <- "not_assignable"
    }
    VDJ$sample_id[i] <- paste(VDJ$sample_id[i], short_fb, sep = "_")
  }

  VDJ$sample_id <- as.factor(VDJ$sample_id)
  VDJ$sample_id <- factor(VDJ$sample_id, levels = levels(VDJ$sample_id))

  if(vdj.only == FALSE){
    vgm[[1]] <- VDJ
  }else{
    return(VDJ)
  }


  #Step 3: Integrate new sample_id information into the seurat object (vgm[[2]]) - Only if integrate.in.gex set to TRUE
  if(integrate.in.gex == TRUE){

    vgm[[2]]$original_sample_id <- rep("",ncol(vgm[[2]]))
    vgm[[2]]$sample_id <- as.character(vgm[[2]]$sample_id) #turn sample_id into character to be able to overwrite the original sample_id's

    #Check if orig_barcode in vgm[[1]] is in the same format as in vgm[[2]]! If it's not change the format of orig_barcode in vgm[[1]] to the same format as in vgm[[2]].
    if(nchar(vgm[[2]]@meta.data[["orig_barcode"]][[1]]) != nchar(vgm[[1]][["orig_barcode"]][[1]])){
      if(nchar(vgm[[2]]@meta.data[["orig_barcode"]][[1]]) < nchar(vgm[[1]][["orig_barcode"]][[1]])){
        nr_characters <- nchar(vgm[[2]]@meta.data[["orig_barcode"]][[1]])
        vgm[[1]]$orig_barcode <- substring(vgm[[1]]$orig_barcode, 1, nr_characters)
      }else{
        nr_characters <- nchar(vgm[[1]][["orig_barcode"]][[1]])
        vgm[[2]]@meta.data$orig_barcode <- substring(vgm[[2]]@meta.data$orig_barcode, 1, nr_characters)
      }
    }


    bc_match_index <- NULL

    for (i in 1:ncol(vgm[[2]])){
      #if you have a match vgm[[1]]$orig_barcode %in% vgm[[2]]$orig_barcode[i] will return one TRUE and the rest FALSE --> length(unique()) will be 2
      if(!vgm[[2]]$orig_barcode[[i]] %in% vgm[[1]]$orig_barcode) { #if no match is found only FALSE is returned for each barcode --> there's only false in the unique() list
        bc_match_index[i] <- NA
        vgm[[2]]$sample_id[[i]] <- paste(vgm[[2]]$sample_id[[i]], "not_assignable", sep = "_") #If we have more cells in the gex data, some cells in vgm[[2]] will not have a barcode match in vgm[[1]]. To those simply assign the original sample_id plus not_assignable. Without this step we would get two groups of not assigned cells for each sample ("s1_not_assignable" and NA)
        vgm[[2]]$original_sample_id[[i]] <- vgm[[2]]$sample_id[[i]]
      }else{
        bc_match_index[i] <- which(vgm[[1]]$orig_barcode == vgm[[2]]$orig_barcode[i]) #which() returns index of the matching vgm[[1]] original barcode
        vgm[[2]]$original_sample_id[[i]] <- as.character(vgm[[1]]$original_sample_id[bc_match_index[i]]) #store original sample_id information
        vgm[[2]]$sample_id[[i]] <- as.character(vgm[[1]]$sample_id[bc_match_index[i]])
      }
    }

    vgm[[2]]$sample_id <- as.factor(vgm[[2]]$sample_id) #Turn sample_id back into a factor
    vgm[[2]]$sample_id <- factor(vgm[[2]]$sample_id, levels = levels(VDJ$sample_id))
  }

  return(vgm)
}
