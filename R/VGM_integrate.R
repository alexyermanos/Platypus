#' (Re)-intergrated VDJ and GEX of one or two separate VGM objects. This can be used as a simple "updating" utility function, if metadata has been added to the VDJ dataframe and is also needed in the GEX matrix or the reverse. Entries are integrated by barcode. If barcodes have been altered (barcode column in VDJ and cell names in GEX), the function will not yield results
#' @param VGM Output object from the VDJ_GEX_matrix function (VDJ_GEX_matrix.output)
#' @param columns.to.transfer Optional. Character Vector. Column names of either the VDJ matrix or GEX meta.data that should be transferred to the corresponding other matrix. if not provided all columns missing from one will be integrated into the other matrix
#' @return An output object from the VDJ_GEX_matrix function with added columns in VDJ or GEX
#' @export
#' @examples
#'
#' #Adding a new clonotyping method to VDJ
#' small_vgm[[1]] <- VDJ_clonotype_v3(VDJ=Platypus::small_vgm[[1]],
#' clone.strategy="cdr3.nt",
#' hierarchical = "single.chains", global.clonotype = TRUE)
#'
#' small_vgm <- VGM_integrate(
#' VGM = small_vgm,
#' columns.to.transfer = NULL) #transfer all new columns
#'
#' small_vgm <- VGM_integrate(
#' VGM = small_vgm,
#' columns.to.transfer = c("global_clonotype_id_cdr3.nt"))
#' #transfer only selected columns

VGM_integrate <- function(VGM,
                          columns.to.transfer){

  if(missing(VGM)) stop("Please provide a VDJ_GEX_matrix output object ")
  if(missing(columns.to.transfer)) columns.to.transfer <- "all"
  if(is.null(columns.to.transfer)) columns.to.transfer <- "all"

  #VERSION is set for now:
  platypus.version <- "v3"

  if(class(VGM[[1]]) != "data.frame"){
    stop("No VDJ matrix found in the provided VGM[[1]]")
  }
  if(!class(VGM[[2]]) %in% c("Seurat","SeuratObject")){
    stop("No GEX object found in the provided VGM[[2]]")
  }

  #checking if there is any overlap between VDJ and GEX
  message(paste0("Cells in GEX: ", ncol(VGM[[2]]),"; cells in VDJ: ", nrow(VGM[[1]]), "\n Overlap: ", length(intersect(colnames(VGM[[2]]), VGM[[1]]$barcode)), " \n"))

  if(length(intersect(colnames(VGM[[2]]), VGM[[1]]$barcode)) == 0) stop("No barcode overlap between VDJ barcode column and GEX colnames. Please verify that barcodes or cell names have not been altered")

  #now find columns to integrate
  #starting with VDJ to GEX
  if(columns.to.transfer != "all"){
    if(any(columns.to.transfer %in% names(VGM[[1]]))){
     cl_VDJ_to_GEX <- columns.to.transfer[columns.to.transfer %in% names(VGM[[1]])]
    } else {
      cl_VDJ_to_GEX <- c()
      }
    } else {
    cl_VDJ_to_GEX <- names(VGM[[1]])[!names(VGM[[1]]) %in% names(VGM[[2]]@meta.data)]
    message("Integrating all VDJ columns which are not in GEX to GEX")
  }
  if(length(cl_VDJ_to_GEX) > 0){
    VGM[[2]]@meta.data[,"int_merge_bc"] <- colnames(VGM[[2]])
    VGM[[1]][,"int_merge_bc"] <- VGM[[1]]$barcode

    to_merge <- VGM[[1]][,c("int_merge_bc",cl_VDJ_to_GEX)]
    new_GEX_meta <- merge(VGM[[2]]@meta.data, to_merge, by = "int_merge_bc", all.x = T, all.y = F, sort = F)
    rownames(new_GEX_meta) <- new_GEX_meta$int_merge_bc
    new_GEX_meta <- new_GEX_meta[,names(new_GEX_meta) != "int_merge_bc"]
    VGM[[2]]@meta.data <- new_GEX_meta
    message(paste0("Added columns to GEX: ", paste0(cl_VDJ_to_GEX, collapse = "; ")))
  }

  #now GEX to VDJ
  if(columns.to.transfer != "all"){
    if(any(columns.to.transfer %in% names(VGM[[2]]@meta.data))){
      cl_GEX_to_VDJ <- columns.to.transfer[columns.to.transfer %in% names(VGM[[2]]@meta.data)]
    } else {
      cl_GEX_to_VDJ <- c()
    }
  } else {
    cl_GEX_to_VDJ <- names(VGM[[2]]@meta.data)[!names(VGM[[2]]@meta.data) %in% names(VGM[[1]])]
    message("Integrating all VDJ columns which are not in GEX to GEX")
  }
  if(length(cl_GEX_to_VDJ) > 0){
    VGM[[2]]@meta.data[,"int_merge_bc"] <- colnames(VGM[[2]])
    VGM[[1]][,"int_merge_bc"] <- VGM[[1]]$barcode

    to_merge <- VGM[[2]]@meta.data[,c("int_merge_bc",cl_GEX_to_VDJ)]
    new_VDJ <- merge(VGM[[1]], to_merge, by = "int_merge_bc", all.x = T, all.y = F, sort = F)
    new_VDJ <- new_VDJ[,names(new_VDJ) != "int_merge_bc"]
    VGM[[1]] <- new_VDJ
    message(paste0("Added columns to VDJ: ", paste0(cl_GEX_to_VDJ, collapse = "; ")))
  }
  return(VGM)
}


