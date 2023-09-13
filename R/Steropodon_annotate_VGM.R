#' Annotates the VGM object with information from the Steropodon object


#'@description Function to integrate information from the Steropodon object (structure clusters, paratope and epitope sequences after docking) into the VGM object (both GEX and VDJ) for downstream analysis.
#' For use cases, see the Steropodon vignette (https://alexyermanos.github.io/Platypus/articles/Steropodon.html) - 9. Cluster analysis via Steropodon_cluster


#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' Steropodon_cluster, Steropodon_dock and other functions should be first ran on the object to obtain the information that will be integrated into the VGM object.
#' @param VGM the Platypus VGM object (GEX and VDJ) to be annotated.

#' @return a VGM object with all information available added from Steropodon objects(s)
#' @export
#' @examples
#' \dontrun{
#' Steropodon_annotate_VGM(steropodon.objkect = structures, VGM = VGM)
#'}



Steropodon_annotate_VGM <- function(steropodon.object,
                                    VGM
                                    ){

  if(missing(steropodon.object)) stop('Please input your Steropodon object/ nested list of objects!')
  if(missing(VGM)) stop('Please input a VGM object or a VDJ/GEX one!')


  annotation_per_barcode <- function(steropodon.object){

    barcodes_df <- as.data.frame(steropodon.object@barcodes)
    colnames(barcodes_df) <- c('barcode')

    if(!is.null(steropodon.object@core)){
      barcodes_df$core <- paste0(bio3d::pdbseq(steropodon.object@core), collapse = '')

    }

    if(!is.null(steropodon.object@cluster)){
      barcodes_df$cluster <- steropodon.object@cluster
    }

    if(!is.null(steropodon.object@epitope)){
      barcodes_df$epitope <- paste0(bio3d::pdbseq(steropodon.object@epitope), collapse = '')

    }

    if(!is.null(steropodon.object@paratope)){
      barcodes_df$paratope <- paste0(bio3d::pdbseq(steropodon.object@paratope), collapse = '')
    }


    return(barcodes_df)
  }

  annotate_VDJ <- function(barcode.df, VDJ){
    cols <- colnames(barcode.df)
    cols <- cols[cols!='barcode']

    for(col in cols){
      VDJ[col] <- NULL
    }

    VDJ <- dplyr::left_join(VDJ, barcode.df, by = c('barcode'))

    return(VDJ)
  }

  annotate_GEX <- function(barcode.df, GEX){
    cols <- colnames(barcode.df)
    cols <- cols[cols!='barcode']

    for(col in cols){
      GEX@meta.data[col] <- NULL
    }

    GEX@meta.data$barcode <- rownames(GEX@meta.data)
    GEX@meta.data <- dplyr::left_join(GEX@meta.data, barcode.df, by = c('barcode'))

    GEX <- Seurat::AddMetaData(object = GEX, metadata = GEX@meta.data)


    return(GEX)
  }

  steropodon_list <- unnest_steropodon(steropodon.object)
  barcode_dfs <- lapply(steropodon_list, function(x) annotation_per_barcode(x))
  barcode_df <- do.call('rbind', barcode_dfs)

  if(inherits(VGM, 'list')){
    VGM[[1]] <- annotate_VDJ(barcode_df, VGM[[1]])
    VGM[[2]] <- annotate_GEX(barcode_df, VGM[[2]])

  }else if(inherits(VGM, 'data.frame')){
    VGM <- annotate_VDJ(barcode_df, VGM)

  }else if(inherits(VGM, 'Seurat')){
    VGM <- annotate_GEX(barcode_df, VGM)

  }else{
    stop('Could not find a valid VGM/VDJ/GEX object to be annotated!')
  }


  return(VGM)
}
