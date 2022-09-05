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
