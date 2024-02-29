#' Utility for VDJ GEX matrix to integrated VDJ and GEX objects after addition of data to either
#'
#'@description (Re)-intergrated VDJ and GEX of one or two separate VGM objects. This can be used as a simple "updating" utility function, if metadata has been added to the VDJ dataframe and is also needed in the GEX matrix or the reverse. Entries are integrated by barcode. If barcodes have been altered (barcode column in VDJ and cell names in GEX), the function will not yield results
#' @param VGM Output object from the VDJ_GEX_matrix function (VDJ_GEX_matrix.output)
#' @param columns.to.transfer Optional. Character Vector. Column names of either the VDJ matrix or GEX meta.data that should be transferred to the corresponding other matrix. if not provided all columns missing from one will be integrated into the other matrix
#' @param genes.to.VDJ Character vector of gene names in GEX. In many cases it is useful to extract expression values for a gene to metadata. This is done via SeuratObject::FetchData(vars  = genes,slot = seurat.slot) function. The VGM integrate takes gene ids, extracts these and adds them to the VDJ dataframe. If provided, no other columns are integrated between VDJ and GEX and columns.to.transfer is ignored.
#' @param seurat.slot GEX object data slot to pull from. Can be 'counts', 'data', or 'scale.data'
#' @return An output object from the VDJ_GEX_matrix function with added columns in VDJ or GEX
#' @export
#' @examples
#'
#' small_vgm[[1]] <- VDJ_clonotype(VDJ=Platypus::small_vgm[[1]],
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
#'
#' small_vgm <- VGM_integrate(
#' small_vgm, genes.to.VDJ = c("CD19","CD24A"),seurat.slot = "counts")

VGM_integrate <- function(VGM,
                          columns.to.transfer,
                          genes.to.VDJ,
                          seurat.slot){

  if(missing(VGM)) stop("Please provide a VDJ_GEX_matrix output object ")
  if(missing(columns.to.transfer)) columns.to.transfer <- "all"
  if(is.null(columns.to.transfer)) columns.to.transfer <- "all"

  #VERSION is set for now:
  platypus.version <- "v3"

  if(!inherits(VGM[[1]],"data.frame")){
    stop("No VDJ matrix found in the provided VGM[[1]]")
  }
  if(!class(VGM[[2]]) %in% c("Seurat","SeuratObject")){
    stop("No GEX object found in the provided VGM[[2]]")
  }

  #checking if there is any overlap between VDJ and GEX
  message(paste0("Cells in GEX: ", ncol(VGM[[2]]),"; cells in VDJ: ", nrow(VGM[[1]]), "\n Overlap: ", length(intersect(colnames(VGM[[2]]), VGM[[1]]$barcode)), " \n"))

  if(length(intersect(colnames(VGM[[2]]), VGM[[1]]$barcode)) == 0) stop("No barcode overlap between VDJ barcode column and GEX colnames. Please verify that barcodes or cell names have not been altered")

  if(missing(genes.to.VDJ) == TRUE){ #if no genes to be extracted to VDJ

  #Keeping original GEX names.
  #As we add VDJ -> GEX first, we need to keep track of original columns in GEX so no "backintegration" happens
  ori_GEX_names <- names(VGM[[2]]@meta.data)
  ori_VDJ_names <- names(VGM[[1]])

  #now find columns to integrate
  #starting with VDJ to GEX
  if(columns.to.transfer[1] != "all"){
    if(any(columns.to.transfer %in% ori_VDJ_names)){
     cl_VDJ_to_GEX <- columns.to.transfer[columns.to.transfer %in% ori_VDJ_names]
     cl_VDJ_to_GEX <- cl_VDJ_to_GEX[!cl_VDJ_to_GEX %in% ori_GEX_names] #no need to transfer columns already in VDJ
    } else {
      cl_VDJ_to_GEX <- c()
      }
    } else {

      VGM[[2]]@meta.data <- VGM[[2]]@meta.data[-c(which(names(VGM[[2]]@meta.data) %in% c("clonotype_id", "clonotype_frequency")))]#remove columns to allow updating these in case a new clonotyping strategy has been added

    cl_VDJ_to_GEX <- ori_VDJ_names[!ori_VDJ_names %in% ori_GEX_names]
    cl_VDJ_to_GEX <- cl_VDJ_to_GEX[!cl_VDJ_to_GEX %in% c("barcode", "GEX_available")]
    cl_VDJ_to_GEX <- c(cl_VDJ_to_GEX, "clonotype_id", "clonotype_frequency") #append columns to allow updating these in case a new clonotyping strategy has been added
    VGM[[2]]@meta.data <- VGM[[2]]@meta.data[-c(which(names(VGM[[2]]@meta.data) %in% c("clonotype_id", "clonotype_frequency")))]#remove columns to allow updating these in case a new clonotyping strategy has been added
    if(length(cl_VDJ_to_GEX)>0)message("Integrating all VDJ columns which are not in GEX to GEX")
    }

  if(length(cl_VDJ_to_GEX) > 0){
    VGM[[2]]@meta.data[,"int_merge_bc"] <- colnames(VGM[[2]])
    VGM[[1]][,"int_merge_bc"] <- VGM[[1]]$barcode

    to_merge <- VGM[[1]][,c("int_merge_bc",cl_VDJ_to_GEX)]
    new_GEX_meta <- merge(VGM[[2]]@meta.data, to_merge, by = "int_merge_bc", all.x = TRUE, all.y = FALSE, sort = FALSE)
    rownames(new_GEX_meta) <- new_GEX_meta$int_merge_bc
    new_GEX_meta <- new_GEX_meta[,names(new_GEX_meta) != "int_merge_bc"]
    VGM[[1]] <- VGM[[1]][,names(VGM[[1]]) != "int_merge_bc"]
    VGM[[2]]@meta.data <- new_GEX_meta
    message(paste0("Added columns to GEX: ", paste0(cl_VDJ_to_GEX, collapse = "; ")))
  }

  #now GEX to VDJ
  if(columns.to.transfer[1] != "all"){

    if(any(columns.to.transfer %in% ori_GEX_names)){
      cl_GEX_to_VDJ <- columns.to.transfer[columns.to.transfer %in% ori_GEX_names]
      cl_GEX_to_VDJ <- cl_GEX_to_VDJ[!cl_GEX_to_VDJ %in% names(VGM[[1]])] #no need to transfer columns already in VDJ
    } else {
      cl_GEX_to_VDJ <- c()
    }
  } else {
    cl_GEX_to_VDJ <- ori_GEX_names[!ori_GEX_names %in% ori_VDJ_names]
    cl_GEX_to_VDJ <- cl_GEX_to_VDJ[!cl_GEX_to_VDJ %in% c("nCount_RNA", "nFeature_RNA", "VDJ_available", "percent.mt", "RNA_snn_res.0.5", "FB_assignment.y")]
    if(length(cl_GEX_to_VDJ)>0) message("Integrating all GEX columns which are not in VDJ to VDJ")
  }
  if(length(cl_GEX_to_VDJ) > 0){
    VGM[[2]]@meta.data[,"int_merge_bc"] <- colnames(VGM[[2]])
    VGM[[1]][,"int_merge_bc"] <- VGM[[1]]$barcode
    to_merge <- VGM[[2]]@meta.data[,c("int_merge_bc",cl_GEX_to_VDJ)]
    new_VDJ <- merge(VGM[[1]], to_merge, by = "int_merge_bc", all.x = TRUE, all.y = FALSE, sort = FALSE)
    VGM[[2]]@meta.data <- VGM[[2]]@meta.data[,names(VGM[[2]]@meta.data) != "int_merge_bc"]
    new_VDJ <- new_VDJ[,names(new_VDJ) != "int_merge_bc"]

    VGM[[1]] <- new_VDJ
    message(paste0("Added columns to VDJ: ", paste0(cl_GEX_to_VDJ, collapse = "; ")))
  }

  } else { #genes.to.VDJ not missing

    #Fetch genes
    fetched <- SeuratObject::FetchData(VGM[[2]], vars = genes.to.VDJ, slot = seurat.slot)
    #Integrate to VDJ via barcode

    fetched[,"int_merge_bc"] <- rownames(fetched)
    VGM[[1]][,"int_merge_bc"] <- VGM[[1]]$barcode

    new_VDJ <- merge(VGM[[1]], fetched, by = "int_merge_bc", all.x = TRUE, all.y = FALSE, sort = FALSE)
    new_VDJ <- new_VDJ[,names(new_VDJ) != "int_merge_bc"]
    VGM[[1]] <- new_VDJ
    message(paste0("Added columns to VDJ: ", paste0(colnames(fetched)[1:ncol(fetched)-1], collapse = "; ")))

  }
  return(VGM)
}
