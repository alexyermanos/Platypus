#' Create a nested list of longitudinal AntibodyForests objects


#'@description Adds the dynamic slots to a nested list of AntibodyForests objects outputted from AntibodyForests function. Also inverts the nested list (per clonotype per sample instead of per sample per clonotype) - for tracking the evolution of a specific clonotype across multiple timepoints (samples).
#' The timepoints order can be specified in the timepoint.order parameter.
#' The new dynamic graphs  contain all the unique nodes across the timepoints, but with edges created only for a single tree of a given timepoint.
#' The new dynamic slots will be used for downstream analyses - AntibodyForests_metrics(graph.type = 'dynamic') and AntibodyForests_track_nodes.
#' Before running this function, ensure your clonotypes are defined the same way across each timepoint before creating your networks using AntibodyForests (e.g., use the VDJ_call_enclone function or VDJ_clonotype with global.clonotype set to TRUE to ensure clonotype 1 is defined the same across each timepoint, otherwise clonotype1 in timepoint/sample 1 might not correspond to the same clonal definition as clonotype1 in timepoint/sample2).

#' @param trees nested list of AntibodyForests objects, as obtained from the AntibodyForests function. Ensure the clonotype definition is consistent across each timepoint before running this function (and before running AntibodyForests to obtain your trees). Also ensure the timepoint ids are present in the sample_id column of your VDJ/VGM[[1]] object.
#' @param graph.type string - 'tree' will use the usual output of the AntibodyForests function (tree graphs), 'heterogeneous' will use the output of the AntibodyForests_heterogeneous function (bipartite networks for both cells and sequences).
#' @param timepoints.order vector of strings - order of the timepoints in the resulting nested list. For example, output[[1]][[1]] denotes the first clonotype, first timepoint/sample, output[[1]][[2]] denotes the first clonotype, second timepoint/sample,


#' @return nested list of AntibodyForests objects for each clonotype and each sample/timepoint. For example, output[[1]][[2]] denotes the AntibodyForests object of the first clonotype, second timepoint.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_dynamics(trees, graph.type = 'tree', timepoint.order = c('s1', 's2', 's3')
#'}



AntibodyForests_dynamics <- function(trees,
                                     graph.type,
                                     timepoints.order){


  if(missing(trees)) stop('Please input a nested list of AntibodyForests objects!')
  if(missing(graph.type)) graph.type <- 'tree'
  if(missing(timepoints.order)) timepoints.order <- names(trees)


  trees <- trees[timepoints.order]

  shared_clonotypes <- names(trees[[1]])
  for(i in 2:length(trees)){
    shared_clonotypes <- intersect(shared_clonotypes, names(trees[[i]]))
  }


  if(length(shared_clonotypes) == 0){
    stop('Could not find a single clonotype shared across each timepoint')
  }


  output_trees <- vector(mode = 'list', length = length(shared_clonotypes))
  names(output_trees) <- shared_clonotypes

  for(i in 1:length(shared_clonotypes)){
    for(j in 1:length(timepoints.order)){
      output_trees[[i]][[j]] <- trees[[timepoints.order[j]]][[shared_clonotypes[i]]]
    }
    names(output_trees[[i]]) <- timepoints.order
  }


  for(i in 1:length(output_trees)){
    new_vertex_df <- vector(mode = 'list', length = length(output_trees[[i]]))
    new_edge_df <- vector(mode = 'list', length = length(output_trees[[i]]))

    shared_sequences <- output_trees[[i]][[1]]@sequences

    for(j in 1:length(output_trees[[i]])){

      if(graph.type == 'tree'){
        g <- output_trees[[i]][[j]]@tree
      }else if(graph.type == 'heterogeneous'){
        g <- output_trees[[i]][[j]]@heterogeneous
      }

      shared_sequences <- intersect(shared_sequences, output_trees[[i]][[j]]@sequences)
      new_vertex_df[[j]] <- igraph::as_data_frame(output_trees[[i]][[j]]@tree, what = 'vertices')
      new_edge_df[[j]] <- igraph::as_data_frame(output_trees[[i]][[j]]@tree, what = 'edges')
    }

    new_vertex_df <- do.call('rbind', new_vertex_df)
    new_edge_df <- do.call('rbind', new_edge_df)

    new_vertex_df$old_label <- new_vertex_df$label
    new_vertex_df$label <- 1:nrow(new_vertex_df)
    new_vertex_df$timepoints <- new_vertex_df$sample_id

    sample_ids <- unique(new_vertex_df$sample_id)
    ids <- 1:length(sample_ids)
    names(ids) <- sample_ids

    new_vertex_df$timepoints_id <- ids[new_vertex_df$timepoints]


    new_vertex_df$shared[new_vertex_df$network_sequences %in% shared_sequences] <- 'yes'
    new_vertex_df$shared[!(new_vertex_df$network_sequences %in% shared_sequences)] <- 'no'


    for(j in 1:length(output_trees[[i]])){

      if(graph.type == 'tree'){
        g <- output_trees[[i]][[j]]@tree
      }else if(graph.type == 'heterogeneous'){
        g <- output_trees[[i]][[j]]@heterogeneous
      }

      directed <- igraph::is_directed(g)
      timepoint <- output_trees[[i]][[j]]@sample_id
      edge_df <- igraph::as_data_frame(g, what = 'edges')

      label_dict <- new_vertex_df$label[new_vertex_df$timepoints == timepoint]
      names(label_dict) <- new_vertex_df$old_label[new_vertex_df$timepoints == timepoint]

      edge_df$from <- label_dict[edge_df$from]
      edge_df$to <- label_dict[edge_df$to]
      new_g <- igraph::graph_from_data_frame(edge_df, directed = directed, vertices = new_vertex_df)

      output_trees[[i]][[j]]@dynamic <- new_g

    }
  }

  return(output_trees)
}
