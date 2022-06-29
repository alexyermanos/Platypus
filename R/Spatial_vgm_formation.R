#' Addition of the spatial information to the vgm matrix, output of VDJ_GEX_matrix function from Platypus.
#' @param vgm Large list, output of VDJ_GEX_matrix function
#' @param tissue_lowres_image_path Path to file containing the image of the tissue in png format
#' @param scalefactors_json_path Path to a file for convertint pixel positions in the original, full-resolution image to pixel positions in the histological image in json format
#' @param tissue_potions_list_path Path to a text file containing a table with rows that correspond to spots in csv format
#' @param cluster_path Path to 10X Genomic clustering file that is not specific for immune cells, in csv format
#' @param matrix_path Path to the filtered feature barcode matrix containing barcode from fixed list of known-good barcode sequences in the h5 format
#' @return Returns a large list containing the VDJ information, the GEX information and the spatial data.
#' @export
#' @examples
#' \dontrun{
#' tissue_lowres_image_path<-list()
#' tissue_lowres_image_path[[1]]<-c("c:/.../tissue_lowres_image.png")
#' 
#' scalefactors_json_path<-list()
#' scalefactors_json_path[[1]]<-c("c:/.../scalefactors_json.json")
#' 
#' tissue_positions_list_path<-list()
#' tissue_positions_list_path[[1]]<-c("c:/.../tissue_positions_list.csv")
#' 
#' cluster_path<-list()
#' cluster_path[[1]]<-c("c:/.../analysis/clustering/graphclust/clusters.csv")
#' 
#' matrix_path<-list()
#' matrix_path[[1]]<-c("c:/.../filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5")
#' 
#' vgm_spatial<-Spatial_vgm_formation(vgm = vgm_without_spatial_data_and_VDJ,
#' tissue_lowres_image_path = tissue_lowres_image_path,tissue_positions_list_path = tissue_positions_list_path, 
#' scalefactors_json_path = scalefactors_json_path, cluster_path = cluster_path, matrix_path = matrix_path)
#'}

Spatial_vgm_formation<-function(vgm,
                                tissue_lowres_image_path,
                                scalefactors_json_path,
                                tissue_positions_list_path,
                                cluster_path,matrix_path){
  
  if(missing(vgm)) stop("Please provide vgm input for this function")
  if(missing(tissue_lowres_image_path)) stop("Please provide tissue_lowres_image_path input for this function")
  if(missing(scalefactors_json_path)) stop("Please provide scalefactors_json_path input for this function")
  if(missing(tissue_positions_list_path)) stop("Please provide tissue_positions_list_path input for this function")
  if(missing(cluster_path)) stop("Please provide cluster_path,matrix_path input for this function")
  if(missing(matrix_path)) stop("Please provide cluster_path,matrix_path input for this function")
  
  platypus.version <- "v3" 
  
  vgm$spatial<-list()
  for (i in 1:length(sample_names)){
    vgm$spatial$image[[i]]<-read.bitmap(tissue_lowres_image_path[[i]])
    vgm$spatial$scalefactor[[i]] <- fromJSON(file = scalefactors_json_path[[i]])
    vgm$spatial$tissue[[i]] <- read.csv(tissue_positions_list_path[[i]], col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
    vgm$spatial$cluster[[i]] <- read.csv(cluster_path[[i]])
    vgm$spatial$matrix[[i]]<-as.data.frame(t(Read10X_h5(matrix_path[[i]])))
  } 
  return(vgm)
}
