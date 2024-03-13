#' Plot lineage trees created with AntibodyForests
#' @description AntibodyForests_plot takes an AntibodyForests-object and outputs a graphs for a specified sample and clonotype.
#' @param AntibodyForests_object list - AntibodyForests-object as obtained from the AntibodyForests() function
#' @param sample string - denotes the sample that contains the clonotype to be plotted.
#' @param clonotype string - denotes the clonotype from which the B cell lineage tree should be plotted.
#' @param node.color string specifying the node feature to be used for coloring the nodes.
#' If the node has multiple values for this node feature, the resulting nodes will be a pie chart
#' If the node.color parameter is NULL (default), the default node color will be light blue for sequence nodes, grey for intermediate/inferred nodes, and orange for germline nodes.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_plot (AntibodyForests_object,
#'                   sample = "S1",
#'                   clonotype = "clonotype1")
#'}

AntibodyForests_plot <- function(AntibodyForests_object,
                              sample,
                              clonotype,
                              node.color){
  
  #Check input and set defaults
  if(missing(AntibodyForests_object)) stop('Please input a AntibodyForests-object, output of AntibodyForests()')
  if(missing(sample)) stop('Please input a sample name from the AntibodyForests-object')
  if(!(sample %in% names(AntibodyForests_object))) stop("Sample name is not in the AntibodyForests-object")
  if(missing(clonotype)) stop('Please input a clonotype name from the AntibodyForests-object')
  if(!(clonotype %in% names(AntibodyForests_object[[sample]]))) stop("Clonotype name is not in the AntibodyForests-object")
  if(missing(node.color)) node.color <- NULL
  if(!(is.null(node.color)) && !(node.color %in% names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[1]]))) stop("node.color is not in the AntibodyForests-object")
  
  
  # Retrieve igraph object from AntibodyForests object
  tree <- AntibodyForests_object[[sample]][[clonotype]][["igraph"]]
  nodes <- AntibodyForests_object[[sample]][[clonotype]][["nodes"]]
  
  # Define node labels, germline becomes "G" and remove "node" for the rest of the node labels
  igraph::V(tree)$label <- ifelse(igraph::V(tree)$name == "germline", "G", 
                                  ifelse(startsWith(igraph::V(tree)$name, "node"), gsub(pattern = "node", replacement = "", igraph::V(tree)$name), igraph::V(tree)$name))
  
  # Arrange the nodes by using 'germline' node as the root and by directing the tree downwards using the 'igraph::layout_as_tree()' function
  layout <- igraph::layout_as_tree(tree, root = "germline")
  
  if (is.null(node.color)){
    # Assign default node colors
    igraph::V(tree)$color <- ifelse(igraph::V(tree)$label == "G", "orange",
                                    ifelse(startsWith(igraph::V(tree)$name, "node"), "lightblue", "grey"))
    
    # Plot tree
    igraph::plot.igraph(tree, layout = layout,
                        vertex.label.dist = 0,
                        vertex.size = 10, 
                        vertex.label.cex = 0.8,
                        edge.arrow.size = 0.1)
  }else{
    #Create a list with the values to color per sequence per node
    color_values <- lapply(nodes,function(x){x[[node.color]]})[names(V(tree))]
    #Give the value "germline" to the germline node
    color_values$germline <- "germline"
    
    #Change NA values to "NA" for plotting purpose
    color_values <- lapply(color_values, function(x){
      index <- which(is.na(x))
      if (length(index) > 0){
        x[[index]] <- "NA"
      }
      return(x)
    })
    
    #Get the unique values to color
    unique_values <- unique(unlist(color_values))
    
    #Specific color scheme for isotypes/VDJ_cgene
    if (node.color %in% c("isotype", "VDJ_cgene")){
      color_list <- dplyr::case_match(unique_values,
                                      "IgA" ~ "red",
                                       "IgH" ~ "purple",
                                       "IgG" ~ "green",
                                       "IgM" ~ "black",
                                       "IgD" ~ "blue",
                                       "IGHA" ~ "red",
                                       "IGHE" ~ "purple",
                                       "IGHG" ~ "green",
                                       "IGHG1" ~ "green4",
                                       "IGHG2B" ~ "green3",
                                       "IGHG2C" ~ "green2",
                                       "IGHG3" ~ "green1",
                                       "IGHM" ~ "black",
                                       "IGHD" ~ "blue",
                                       "germline" ~ "orange",
                                      "NA" ~ "grey")
    }else{
      color_list <- RColorBrewer::brewer.pal(length(unique_values), "Set1")
    }
    
      #Create list with the number of sequences per feature per node
      color_counts <- lapply(color_values, function(x){
        vector <- c()
        for (index in 1:length(unique_values)){
          count <- sum(str_count(x, pattern = unique_values[index]))
          vector[index] <- count
        }
        return(vector)
      })
    
    
    # Plot tree
    igraph::plot.igraph(tree, layout = layout,
                        vertex.shape = "pie",
                        vertex.pie=color_counts,
                        vertex.pie.color=list(color_list),
                        vertex.label = NA,
                        vertex.size = 10, 
                        vertex.label.cex = 0.8,
                        edge.arrow.size = 0.1)
    # Add a legend
    graphics::legend("topleft", legend = unique_values, fill = color_list, title = node.color,
                     cex = 0.7)      
        
  }

    

  
  

  

}