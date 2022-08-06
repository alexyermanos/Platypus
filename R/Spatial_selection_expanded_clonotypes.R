#' Selection of VGM[[1]]/VDJ data of the x more expanded clonotypes.
#' @param nb_clonotype Number that describes how many clonotypes we want to extract from the VGM[[1]].
#' @param vgm_VDJ Data frame containing VDJ information, found in the vgm made by platypus. It must have x and y coordinates column and the column containing the factor to plot.
#' @return Returns a data frame with only the data belonging to the number of selected clonotypes. The clonotypes being the most expanded ones.
#' @export
#' @examples
#' \dontrun{
#' top_5_VDJ_data<-Spatial_selection_expanded_clonotypes(nb_clonotype = 5, vgm_VDJ = vgm$VDJ)
#' }

Spatial_selection_expanded_clonotypes<-function(nb_clonotype,
                                                     vgm_VDJ){

  if(missing(nb_clonotype)) stop("Please provide nb_clonotype input for this function")
  if(missing(vgm_VDJ)) stop("Please provide vgm_VDJ input for this function")

  platypus.version <- "v3"

  #clonotype
  clonotype<-vgm_VDJ[,c("clonotype_id_10x","clonotype_frequency")]
  names(clonotype)<-c("clonotype_id","frequency")
  clonotype<-clonotype[!duplicated(clonotype), ]
  #selection of cell
  cell<-vgm_VDJ
  clonotype<-clonotype[order(clonotype$frequency, decreasing = T), ]
  b<-list()
  for (i in 1:nb_clonotype) {
    a<-dplyr::filter(vgm_VDJ, vgm_VDJ$clonotype_id==clonotype$clonotype_id[[i]])
    b<-rbind(b,a)
  }
  b<-as.data.frame(b)
  return(b)
}
