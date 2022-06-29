#' Selection in a data frame of the most expanded clonotypes data obtained by a simulation of the VDJ by Echidna.
#' @param nb_clonotype Number that describes how many clonotypes we want to select.
#' @param vgm_VDJ Data frame containing all the data on the cell. It must contain the column clonotype_id which describes the number of the clonotype to which the cell belongs. This data frame can be obtained by the assignment functions (VDJ_assignment_random_based, VDJ_assignment_density_based and VDJ_assignment_germline_based).
#' @param clonotype Data frame containing the simulation of the cells wanted by Echidna. It must contain the clonotype_id column. It can be found in the output named "clonotypes" of Echidna's simulate_repertoire function.
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
    a<-filter(vgm_VDJ, vgm_VDJ$clonotype_id==clonotype$clonotype_id[[i]])
    b<-rbind(b,a)
  }
  return(b)
}
