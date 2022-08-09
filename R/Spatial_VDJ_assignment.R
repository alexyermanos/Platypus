#' Assign simulated immune repertoire sequences (BCR or TCR) simulated by Echidna to transcriptome and location in a spatial image in function of cell type. 
#' @param GEX_matrix Dataframe containing barcode, imagecol and imagerow from bcs_merge. 
#' @param vgm Output of VDJ_GEX_matrix function with already the simulated VDJ data.
#' @param vgm_VDJ Dataframe from VDJ_GEX_matrix output (vgm[[1]]).
#' @param celltype Character designating the cell type that we want to study either "B" or "T".
#' @param simulated_VDJ Large list, output of Echidna simulate_repertoire function. Only needed if we use simulated data.
#' @param method Character to chose the assignment method of BCR or TCR to transcriptomic information, it can be "random", "density" or "germline". 
#' @return A dataframe corresponding to the VDJ (VGM[[1]]) with GEX_barcode and x, y coordinates column (allowing to localise each BCR or TCR on the spatial image).
#' @export
#' @examples
#' \dontrun{
#' #1)Assignment random to GEX
#' random_BCR_assignment <- Spatial_VDJ_assignment(GEX_matrix = GEX_matrix,vgm = vgm_with_simulated_VDJ,
#' vgm_VDJ = vgm_with_simulated_VDJ$VDJ, celltype = "B", simulated_VDJ = simulated_B_cells_VDJ, method = "random")
#' 
#' #2)Assignment density-based
#' density_BCR_assignment<-Spatial_VDJ_assignment(GEX_matrix = GEX_matrix,vgm = vgm_with_simulated_VDJ,
#' vgm_VDJ = vgm_with_simulated_VDJ$VDJ, celltype = "B", simulated_VDJ = simulated_B_cells_VDJ, method = "density")
#' 
#' #3)Assignment germline-based
#' germline_BCR_assignment<-Spatial_VDJ_assignment(GEX_matrix = GEX_matrix,vgm = vgm_with_simulated_VDJ,
#' vgm_VDJ = vgm_with_simulated_VDJ$VDJ, celltype = "B", simulated_VDJ = simulated_B_cells_VDJ, method = "germline")
#'}

Spatial_VDJ_assignment<-function(GEX_matrix, 
                                 vgm,
                                 vgm_VDJ,
                                 celltype,
                                 simulated_VDJ, 
                                 method = c("random","density","germline")){
  
  if(missing(GEX_matrix)) stop("Please provide GEX_matrix input for this function")
  if(missing(vgm)) stop("Please provide vgm input for this function")
  if(missing(celltype)) stop("Please provide celltype input for this function")
  if(missing(simulated_VDJ)) stop("Please provide simulated_VDJ input for this function")
  if(missing(method)) stop("Please provide method input for this function: random, density or germline")
  
  platypus.version <- "v3" 
  
  raw_clonotype_id = NULL
  barcode = NULL
  clonotype_id_10x = NULL
  clonotype_frequency = NULL
  GEX_barcode = NULL
  frequency = NULL
  
  
  filtered_contig_annotation<-simulated_VDJ$all_contig_annotations
  clonotype<-simulated_VDJ$clonotypes
  #GEX matrix formation containing GEX barcode, imagecol and imagerow
  cell.state<-list()
  cell.state$cell.state<-vgm$GEX$cell.state
  cell.state$barcode<-vgm$GEX$orig_barcode
  cell.state<-as.data.frame(cell.state)
  GEX_matrix$barcode<-gsub("-1","",as.character(cell.state$barcode))
  GEX_matrix<-merge(GEX_matrix,cell.state, by="barcode")
  filtered_rows_cells <- GEX_matrix$cell.state == celltype
  GEX_matrix <- GEX_matrix[filtered_rows_cells, ]
  names(GEX_matrix)[1]<-"GEX_barcode"
  names(GEX_matrix)[2]<-"x"
  names(GEX_matrix)[3]<-"y"
  #Germline assignment, the first cell of each clonotype in the VDJ simulation is determined as the germline cell
  germline<-dplyr::select(filtered_contig_annotation,raw_clonotype_id, barcode)
  nb_clonotypes<-2*length(clonotype$clonotype_id)
  germline_or_clone<-list()
  for(i in 1:length(germline$raw_clonotype_id)){
    if (i <= nb_clonotypes){
      germline_or_clone[[i]]<-"germline"
    }else{
      germline_or_clone[[i]]<-"clone"
    }
  }
  germline_or_clone<-as.data.frame(germline_or_clone)
  germline_or_clone<-t(germline_or_clone)
  germline_or_clone<-as.data.frame(germline_or_clone)
  germline_or_clone<-dplyr::bind_cols(germline,germline_or_clone)
  germline_or_clone$barcode<-gsub("-1","",as.character(germline_or_clone$barcode))
  names(germline_or_clone)[2]<-"orig_barcode"
  #VDJ matrix, VDJ barcode and genetic information
  VDJ_simulated<-vgm_VDJ
  VDJ_simulated<-merge(VDJ_simulated,germline_or_clone, by = "orig_barcode")
  VDJ_simulated<-VDJ_simulated[!duplicated(VDJ_simulated$orig_barcode),]
  names(VDJ_simulated)[50]<-"germline_or_clone"
  VDJ_matrix<-dplyr::select(VDJ_simulated, clonotype_id_10x, clonotype_frequency, barcode, germline_or_clone)
  names(VDJ_matrix)[1]<-"clonotype"
  names(VDJ_matrix)[2]<-"frequency"
  clonotype<-clonotype[order(clonotype$frequency, decreasing = T), ]
  
  #Different method of assignment
  if (method == "random"){
    available_cells<-GEX_matrix
    available_VDJ<-VDJ_matrix
    assignment_matrix<-dplyr::bind_cols(GEX_matrix,VDJ_matrix)
    assignment_matrix<-dplyr::select(assignment_matrix,GEX_barcode,x,y, frequency, barcode)
    VDJ_simulated<-merge(VDJ_simulated, assignment_matrix,by="barcode")
    results<-VDJ_simulated
  }
  if (method == "density"){
    available_cells<-GEX_matrix
    available_VDJ<-VDJ_matrix
    assignment_matrix<-list()
    a<-list()
    for (i in 1:length(clonotype$clonotype_id)){
      if (length(available_cells$x)!=1){
        nb_clonotypes_x<-clonotype$frequency[[i]]
        proba_density_x<-stats::density(available_cells$x)
        proba_x<- stats::approxfun(proba_density_x$x, proba_density_x$y) #approximation function for density function
        proba_density_x<-proba_x(available_cells$x) #find density value from density function for each x coordinates
        proba_density_x<-as.data.frame(proba_density_x)
        #for y
        proba_density_y<-stats::density(available_cells$y)
        proba_y<- stats::approxfun(proba_density_y$x, proba_density_y$y) #approximation function for density function
        proba_density_y<-proba_y(available_cells$y) #find density value from density function for each y coordinates
        proba_density_y<-as.data.frame(proba_density_y)
        #add to available_T_cells
        used_cells<-dplyr::bind_cols(available_cells,proba_density_x,proba_density_y)
        #calcul all distances between each cell and germline cell
        used_cells<-used_cells[order(used_cells$proba_density_x & used_cells$proba_density_y, decreasing = T), ]
        #assignment of VDJ barcode to clone cells
        subset_used_cells<-used_cells[1:nb_clonotypes_x,]
        subset_available_VDJ<-available_VDJ[available_VDJ$clonotype==clonotype$clonotype_id[[i]],]
        subset_available_VDJ<-subset_available_VDJ[order(subset_available_VDJ$germline_or_clone, decreasing = T), ]
        a<-dplyr::bind_cols(subset_used_cells,subset_available_VDJ)
        assignment_matrix<-dplyr::bind_rows(assignment_matrix,a)
        #remove cell already used
        for(n in 1:length(subset_used_cells$GEX_barcode)){
          available_cells<-available_cells[!(available_cells$GEX_barcode==subset_used_cells$GEX_barcode[[n]]),]
        }
        a<-list()
        used_cells<-list()
      } else {
        subset_available_VDJ<-available_VDJ[available_VDJ$clonotype==clonotype$clonotype_id[[i]],]
        a<-dplyr::bind_cols(available_cells,subset_available_VDJ)
        assignment_matrix<-dplyr::bind_rows(assignment_matrix,a)
        available_cells<-available_cells[-(available_cells$GEX_barcode==subset_used_cells$GEX_barcode[[1]]),]
      }
    }
    #assignment_T_matrix 
    assignment_matrix<-dplyr::select(assignment_matrix,GEX_barcode,x,y, frequency, barcode)
    row.names(assignment_matrix)<-c(1:length(assignment_matrix$x))
    #merge with VDJ_simulated
    VDJ_simulated<-merge(VDJ_simulated,assignment_matrix,by="barcode")
    results<-VDJ_simulated
  }
  if (method == "germline"){
    available_cells<-GEX_matrix
    available_VDJ<-VDJ_matrix
    assignment_matrix<-list()
    a<-list()
    for (i in 1:length(clonotype$clonotype_id)){
      if (length(available_cells$x)!=1){
        nb_clonotypes_x<-clonotype$frequency[[i]]
        proba_density_x<-stats::density(available_cells$x)
        proba_x<- stats::approxfun(proba_density_x$x, proba_density_x$y) #approximation function for density function
        proba_density_x<-proba_x(available_cells$x) #find density value from density function for each x coordinates
        proba_density_x<-as.data.frame(proba_density_x)
        #for y
        proba_density_y<-stats::density(available_cells$y)
        proba_y<- stats::approxfun(proba_density_y$x, proba_density_y$y) #approximation function for density function
        proba_density_y<-proba_y(available_cells$y) #find density value from density function for each y coordinates
        proba_density_y<-as.data.frame(proba_density_y)
        #add to available_T_cells
        used_cells<-dplyr::bind_cols(available_cells,proba_density_x,proba_density_y)
        #Order available_T_cells with highest proba density
        #used_cells<-used_cellss[order(used_cells$proba_density_x & used_cells$proba_density_y, decreasing = T), ]
        #cell chosen as germline
        x<-used_cells[which.max(used_cells$proba_density_x), ]
        x_max<-x$x
        y<- used_cells[which.max(used_cells$proba_density_y), ]
        y_max<-y$y
        #calcul all distances between each cell and germline cell
        for (j in 1:length(used_cells$x)){
          used_cells$x1[[j]]<-abs(used_cells$x[[j]]-x_max)
          used_cells$y1[[j]]<-abs(used_cells$y[[j]]-y_max)
          used_cells$distance[[j]]<-sqrt((used_cells$x1[[j]])^2+(used_cells$y1[[j]])^2)
        }
        used_cells$distance<-as.numeric(used_cells$distance)
        used_cells<-used_cells[order(used_cells$distance, decreasing = F), ] #class cell from closest to farthest
        #assignment of VDJ barcode to clone cells
        subset_used_cells<-used_cells[1:nb_clonotypes_x,]
        subset_available_VDJ<-available_VDJ[available_VDJ$clonotype==clonotype$clonotype_id[[i]],]
        subset_available_VDJ<-subset_available_VDJ[order(subset_available_VDJ$germline_or_clone, decreasing = T), ]
        a<-dplyr::bind_cols(subset_used_cells,subset_available_VDJ)
        assignment_matrix<-dplyr::bind_rows(assignment_matrix,a)
        #remove cell already used
        for(n in 1:length(subset_used_cells$GEX_barcode)){
          available_cells<-available_cells[!(available_cells$GEX_barcode==subset_used_cells$GEX_barcode[[n]]),]
        }
        a<-list()
        used_cells<-list()
      } else {
        subset_available_VDJ<-available_VDJ[available_VDJ$clonotype==clonotype$clonotype_id[[i]],]
        a<-dplyr::bind_cols(available_cells,subset_available_VDJ)
        assignment_matrix<-dplyr::bind_rows(assignment_matrix,a)
        available_cells<-available_cells[-(available_cells$GEX_barcode==subset_used_cells$GEX_barcode[[1]]),]
      }
    }
    #assignment_T_matrix 
    assignment_matrix<-dplyr::select(assignment_matrix,GEX_barcode,x,y, frequency, barcode)
    row.names(assignment_matrix)<-c(1:length(assignment_matrix$x))
    #merge with VDJ_simulated
    VDJ_simulated<-merge(VDJ_simulated,assignment_matrix,by="barcode")
    results<-VDJ_simulated
  }
  return(results)
}
