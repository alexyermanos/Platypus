#' This function plots pseudotime along the trajectories which have been constructed with the GEX_trajectories() function.
#' @param cds cell data set object. Output element [[1]] of the GEX_trajectories() function
#' @param root.nodes For monocle3: Root nodes to determine for the pseudotime trajectories. GEX_trajectories output [[3]] yields all the possible root nodes. Choose the ones you like. 
#' @param root.state For monocle2: Root state to determine starting cluster for the pseudotime trajectories. GEX_trajectories ouput [[3]] yields all the possible root states Choose the one you like.
#' @param monocle.version Version of monocle. Either monocle2 or monocle3. Has to be the same as in GEX_trajectories().Default is monocle3. 
#' @return Returns a list.Element [[1]] cell data set object with the pseudotime trajectories. Element [[2]] pseudotime trajectory plot
#' @export
#' @examples
#' \dontrun{ 
#' ##monocle3
#' pseudotime_output <- GEX_pseudotime_trajectory_plot(GEX_trajectories_output[[1]], root.nodes = c('Y_742','Y_448','Y_964') )
#' 
#' ##monocle2
#' pseudotime_output <- GEX_pseudotime_trajectory_plot(GEX_trajectories_output[[1]], monocle.version = 'monocle3', root.state = "2")
#'  }

GEX_pseudotime_trajectory_plot <- function(cds,
                                           root.nodes,
                                           monocle.version,
                                           root.state){
  
  if(missing(monocle.version)){ monocle.version <- "monocle3"}
  if(missing(cds)) stop("Please provide cds object for this function. Output element [[1]] from the GEX_trajectories() function.")
  
  platypus.version <- "v3"
  
  if (monocle.version == "monocle3"){
    if(missing(root.nodes)) stop("Please provide a vector or a single root node, infered out of the GEX_trajectories() function")
    
    print("calculating pseudotime trajectories for each cell depending on the chosen root nodes")
    cds <- monocle3::order_cells(cds = cds, root_pr_nodes = root.nodes) #calculates pseudotime for each cell 

    
    pseudotime.plot <- monocle3::plot_cells(cds,
                                      color_cells_by = "pseudotime",
                                      label_cell_groups=FALSE,
                                      label_leaves=FALSE,
                                      label_branch_points=FALSE,
                                      graph_label_size=1.5) #plot cells colored in pseudotime according to their ordering
    
    #add pseudotime column to cds to later update VGM[[2]]
    SummarizedExperiment::colData(cds)$pseudotime <- monocle3::pseudotime(cds, reduction_method = "UMAP")
    
    return(list(cds, pseudotime.plot))
    
    
  }else if(monocle.version == "monocle2"){
    if(missing(root.state)){ root.state <- "1" }
    
    pseudotime.plot <- monocle::plot_cell_trajectory(cds, color_by = "Pseudotime", root_state = root.state)
  
    #monocle2 uses CellDataSet object (monocle3 uses cell data set object). Here, pseudotime is already as a column integrated in cds
    
    return(list(cds, pseudotime.plot))
  }
  
}




