#' Adds discrete columns containing TRUE / FALSE on whether a given cell is part of a expanded or not-expanded clonotype. Threshold frequency can be set.
#' @param VGM Output object from the VDJ_GEX_matrix function (VDJ_GEX_matrix.output)
#' @param add.to.VDJ Boolean. Whether to add expanded columns to VDJ matrix. Defaults to TRUE
#' @param add.to.GEX Boolean. Whether to add expanded columns to GEX matrix. Defaults to TRUE
#' @param expansion.threshold Integer. Defaults to 1. Cells in clonotypes above this threshold will be marked as expanded = TRUE.
#' @return An output object from the VDJ_GEX_matrix function with added columns containing TRUE / FALSE values based on clonotype frequency.
#' @export
#' @examples
#'
#' #Add info to whole VGM object
#' VGM <- VGM_expanded_clones(
#' VGM = Platypus::small_VGM, add.to.VDJ = TRUE, add.to.GEX = TRUE,
#' expansion.threshold = 1)
#'

VGM_expanded_clones <- function(VGM,
                                add.to.VDJ,
                                add.to.GEX,
                                expansion.threshold){

  if(missing(add.to.VDJ)) add.to.VDJ <- T
  if(missing(add.to.GEX)) add.to.GEX <- T
  if(missing(expansion.threshold)) expansion.threshold <- T

  #VERSION is set for now:
  platypus.version <- "v3"

  if(class(VGM) != "list"){
    stop("For VGM please input a complete VDJ_GEX_matrix object (list class)")
  }
  if(add.to.VDJ & class(VGM[[1]]) == "character"){
    warning("Expanded column could not be added to VDJ matrix, as no VDJ matrix was found in this VGM object")
  } else if(add.to.VDJ & class(VGM[[1]]) != "character"){
    #add expanded columns

    #search clonotype frequency columns (this is in conjunction with the clonotyping V3 function that can add multiple clonotyping strategies to the same VGM)
    freq_cols <- which(stringr::str_detect(names(VGM[[1]]), "clonotype_frequency"))
    new_cols <- c()
    for(i in freq_cols){
      if(class(VGM[[1]][,i]) %in% c("numeric", "integer")){
        VGM[[1]][,paste0("expanded_", names(VGM[[1]])[i])] <- VGM[[1]][,i] > expansion.threshold
        new_cols <- c(new_cols, paste0("expanded_", names(VGM[[1]])[i]))
      } else {
        warning(paste0("Column ",names(VGM[[1]])[i], " is not numeric and will be skipped"))
      }
    }
    message(paste0("Added new columns to VDJ: ", paste0(new_cols, collapse = "; ")))
  }


  if(add.to.GEX & class(VGM[[2]]) == "character"){
    warning("Expanded column could not be added to GEX matrix, as no GEX matrix was found in this VGM object")
  } else if(add.to.GEX & class(VGM[[2]]) != "character"){
    #add expanded columns

    #search clonotype frequency columns (this is in conjunction with the clonotyping V3 function that can add multiple clonotyping strategies to the same VGM)
    freq_cols <- which(stringr::str_detect(names(VGM[[2]]@meta.data), "clonotype_frequency"))
    new_cols <- c()
    for(i in freq_cols){
      if(class(VGM[[2]]@meta.data[,i])  %in% c("numeric", "integer")){
        VGM[[2]]@meta.data[,paste0("expanded_", names(VGM[[2]]@meta.data)[i])] <- VGM[[2]]@meta.data[,i] > expansion.threshold
        new_cols <- c(new_cols, paste0("expanded_", names(VGM[[2]]@meta.data)[i]))
      } else {
        warning(paste0("Column ",names(VGM[[2]]@meta.data)[i], " is not numeric and will be skipped"))
      }
    }
    message(paste0("Added new columns to GEX: ", paste0(new_cols, collapse = "; ")))
  }
  return(VGM)
}


