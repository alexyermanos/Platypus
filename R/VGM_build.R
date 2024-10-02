#' Function to obtain the VGM object by integrating the VDJ and GEX/Seurat objects
#'@description Takes as input a VDJ data frame (as obtained from the VDJ function in Platypus) and a Seurat object. Outputs an integrated VGM object (a list with the first element - the VDJ object; second element - the Seurat object).
#' Integration involves matching by barcodes and adding all features from one object to the other and vice-versa.
#' Authors: Tudor-Stefan Cotet, Victor Kreiner
#' @param VDJ VDJ data frame, obtained from the Platypus VDJ() function
#' @param GEX Seurat object for the single-cell immune receptor repertoire analysis.
#' @param merge.by string - the column name to match both objects/dataframes by. Should be present in both objects (Seurat object meta.data and VDJ dataframe). Defaults to "barcode."
#' @param additional.dataframes vector of data frames - additional dataframes to be matched/merged to the VDJ and GEX. Will be matched by the column denoted in the merge.by parameter (should be present in the VDJ and all subsequent dataframes).
#' @param columns.to.transfer string or vector of strings - columns that should be transferred/appended across all objects (matched by the merge.by parameter). Defaults to "all" - all unique columns from GEX not present in VDJ and vice-versa.
#' @return An output VGM object: a list with the first element - the VDJ object; second element - the GEX/Seurat object. Additonal elements are appended to the list if additional.dataframes is not null.
#' @export
#' @examples
#'\donttest{
#'try({
#' small_vgm <- VGM_build(
#' VDJ = small_vgm[[1]],
#' GEX = small_vgm[[2]],
#' columns.to.transfer = 'all') #transfer all new columns
#' })
#'}


VGM_build <- function(VDJ,
                GEX,
                merge.by,
                additional.dataframes,
                columns.to.transfer
                ){

    if(missing(VDJ)) stop('Please input your VDJ data frame, obtained from the VDJ() function')
    if(missing(GEX)) GEX <- NULL
    if(missing(additional.dataframes) & !is.null(GEX)) additional.dataframes <- NULL
    if(missing(additional.dataframes) & is.null(GEX)) stop('Please input your data frames that need to be merged to the VDJ object')
    if(missing(merge.by)) merge.by <- 'barcode'
    if(missing(columns.to.transfer)) columns.to.transfer <- 'all'


    extract_columns <- function(from.df, to.df, columns.to.transfer){
      # Subroutine to extract unique columns to be added to the to.df
      # Arguments:
      # - from.df: data frame, from which features (columns.to.transfer) will be transferred.
      # - to.df: data frame, target data frame to which the new columns will be transferred.
      # - columns.to.transfer: vector or "all", columns to be transferred. If "all", will append all unique columns in the from.df that are not present in the to.df.
      # Output: vector of column names to be transferred.
      # Authors: Tudor-Stefan Cotet, Victor Kreiner

      # If "all", will add all unique columns in from.df that are not present in to.df. Else will transfer the column names specified in columns.to.transfer
      if(columns.to.transfer[1] != "all"){

        if(any(columns.to.transfer %in% names(from.df))){
          # Ensure the columns in columns.to.transfer are present in the colnames of form.df
          from_to_columns <- columns.to.transfer[columns.to.transfer %in% names(from.df)]
          # No need to transfer columns already in to.df
          from_to_columns <- from_to_columns[!from_to_columns %in% names(to.df)]
        }else{
          # Will return an empty vector if no valid column names have been found.
          from_to_columns <- c()
        }
      }else{
        # Pick unique column names from from.df not present in to.df if columns.to.transfer is set to "all".
        from_to_columns <- names(from.df)[!names(from.df) %in% names(to.df)]
      }

      return(from_to_columns)
    }

    simple_merge <- function(from.df, to.df, by, columns){
      # Simple merging subroutine, matching columns from from.df to to.df by IDs in a given column.
      # Arguments:
      # - from.df: data frame, from which features (specified columns) will be transferred.
      # - to.df: data frame, target data frame to which the new columns will be transferred.
      # - by: string, columns by which merging will performed.
      # - columns: vector of column names to be transferred, as obtained from the extract_columns function.
      # Output: merged data frame with the new columns name.
      # Authors: Tudor-Stefan Cotet

      # Will only merge if the columns parameter is not empty
      if(length(columns) > 0){
        # Subset the data frame that will be added/merged by only the columns to merge and column that will be the merging pivot
        to_merge <- from.df[, c(by, columns)]
        # Perform merge using R's merge() function. Will append all columns in the to.df and no columns from the from.df
        new_df <- merge(to.df, to_merge, by = by, all.x = TRUE, all.y = FALSE, sort = FALSE)
      }else{
        # If there are no columns to merge, will output the to.df instead (unmerged).
        new_df <- to.df
      }

      return(new_df)
    }

    # Create a list for the VGM object - will include the VDJ object (first item), Seurat object (second item), and subsequent data frames (e.g., per sequence protein language model embeddings as obtained from PlatypusPython)
    VGM <- vector("list", 2 + length(additional.dataframes))
    # Merging additional dataframes to the VDJ object
    if(!is.null(additional.dataframes)){
      # Iterate through all data frames to be merged from the additional.dataframes list
      for(i in 1:length(additional.dataframes)){
        # First merge the VDJ information to the additional data frame
        columns_VDJ_to_df <- extract_columns(VDJ, additional.dataframes[i])
        VGM[[i+2]] <- simple_merge(VDJ, additional.dataframes[i], merge.by, columns_VDJ_to_df)

        # First merge the additional data frame information to the VDJ
        columns_df_to_VDJ <- extract_columns(additional.dataframes[i], VDJ)
        VGM[[1]] <- simple_merge(additional.dataframes[i], VDJ, merge.by, columns_VDJ_to_df)
      }
    }

    # If the VGM list already has info, will select the first element as the VDJ. Else VDJ is as specified in the default input argument of the VGM() function
    if(length(VGM[[1]]) > 0){
      VDJ <- VGM[[1]]
    }

    # Merging the VDJ to the GEX object and vice-versa. Adapted from the old VGM_integrate function.
    if(!is.null(GEX)){
      # Extract Seurat meta.data data frame for merging
      GEX_df <- GEX@meta.data

      # Append barcodes to meta.data Seurat df in case of merging by barcodes
      GEX_df$barcode <- rownames(GEX_df)
      GEX_df$barcode <- gsub(GEX_df$barcode, pattern="-1", replacement="")

      # Merge VDJ to GEX
      cl_VDJ_to_GEX <- extract_columns(VDJ, GEX_df, columns.to.transfer)
      # These column names are redundant/ not needed
      cl_VDJ_to_GEX <- cl_VDJ_to_GEX[!cl_VDJ_to_GEX %in% c("barcode", "GEX_available")]

      # Append columns to allow updating these in case a new clonotyping strategy has been added
      cl_VDJ_to_GEX <- c(cl_VDJ_to_GEX, "clonotype_id", "clonotype_frequency")
      GEX_df <- GEX_df[, !names(GEX_df) %in% c("clonotype_id", "clonotype_frequency")]

      # Check if there are columns to merges
      if(length(cl_VDJ_to_GEX) > 0){
        GEX_df <- simple_merge(VDJ, GEX_df, merge.by, cl_VDJ_to_GEX)
      }

      # Merge GEX to VDJ
      # Extract columns in GEX not already present in VDJ
      cl_GEX_to_VDJ <- extract_columns(GEX_df, VDJ, columns.to.transfer)
      # These columns are often redundant
      cl_GEX_to_VDJ <- cl_GEX_to_VDJ[!cl_GEX_to_VDJ %in% c("nCount_RNA", "nFeature_RNA", "VDJ_available", "percent.mt", "RNA_snn_res.0.5", "FB_assignment.y")]

      # Check if there are columns to merge
      if(length(cl_GEX_to_VDJ) > 0){
        VDJ <- simple_merge(GEX_df, VDJ, merge.by, cl_GEX_to_VDJ)
      }

      #Remove duplicate columns (merge still adds them...)
      GEX_df$clonotype_id.1 <- NULL
      GEX_df$clonotype_frequency.1 <- NULL
      VDJ$clonotype_id.1 <- NULL
      VDJ$clonotype_frequency.1 <- NULL

      GEX@meta.data <- GEX_df

    }else{
      message("Could not find GEX/ Seurat object")
    }

  # Add the merged VDJ and GEX to the final VGM object
  VGM[[1]] <- VDJ
  VGM[[2]] <- GEX

  return(VGM)
}
