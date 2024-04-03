#' Plot lineage trees created with AntibodyForests
#' @description AntibodyForests_plot takes an AntibodyForests-object and outputs a graphs for a specified sample and clonotype.
#' @param AntibodyForests_object list - AntibodyForests-object as obtained from the AntibodyForests() function
#' @param sample string - denotes the sample that contains the clonotype to be plotted.
#' @param clonotype string - denotes the clonotype from which the B cell lineage tree should be plotted.
#' @param node.color string specifying the node feature to be used for coloring the nodes.
#' If the node has multiple values for this node feature, the resulting nodes will be a pie chart
#' If the node.color parameter is NULL (default), the default node color will be light blue for sequence nodes, grey for intermediate/inferred nodes, and orange for germline nodes.
#' @param custom.colors string - colors to use for node.color. If the number of colors in custom.colors is less than the number of unique node.color values in a tree, the colors will be repeated.
#' If the custom.color parameter is NULL (default), Set1 of the RColorBrewer package will be used.
#' @param node.size string - If set to "expansion" nodes size relates to the number of cells per node/sequence in this clonotype. If NULL, all nodes have similar size (default).
#' @param node.size.scale integer - scaling factor for the nodes size, default is 10.
#' @param node.label boolean - If TRUE plots then node number. Default is FALSE.
#' @param edge.length string - If set to "distance", the distance calculated by AntibodyForests() will be used as edge length. If NULL (default), all edges will have the same length.
#' @param legend.position string - "topleft" (default), "topright", "bottomleft", or "bottomright"
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
                              node.color,
                              custom.colors,
                              node.size,
                              node.size.scale,
                              node.label,
                              edge.length,
                              legend.position){
  
  #Check input and set defaults
  if(missing(AntibodyForests_object)) stop('Please input a AntibodyForests-object, output of AntibodyForests()')
  if(missing(sample)) stop('Please input a sample name from the AntibodyForests-object')
  if(!(sample %in% names(AntibodyForests_object))) stop("Sample name is not in the AntibodyForests-object")
  if(missing(clonotype)) stop('Please input a clonotype name from the AntibodyForests-object')
  if(!(clonotype %in% names(AntibodyForests_object[[sample]]))) stop("Clonotype name is not in the AntibodyForests-object")
  if(missing(node.color)) node.color <- NULL
  if(!(is.null(node.color)) && !(node.color %in% names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[1]]))) stop("node.color is not in the AntibodyForests-object")
  if(missing(custom.colors)) custom.colors <- NULL
  if(missing(node.size)) node.size <- NULL
  if(missing(node.size.scale)) node.size.scale <- 10
  if(missing(node.label)) node.label <- F
  if(missing(edge.length)) edge.length <- NULL
  if(missing(legend.position)) legend.position <- "topleft"
  
  #1. Retrieve igraph object and node information from AntibodyForests-object
  tree <- AntibodyForests_object[[sample]][[clonotype]][["igraph"]]
  nodes <- AntibodyForests_object[[sample]][[clonotype]][["nodes"]]
  
  #2. Set node labels
  if (node.label){
    #germline becomes "G" and remove "node" for the rest of the node labels
    igraph::V(tree)$label <- ifelse(igraph::V(tree)$name == "germline", "G", 
                                    ifelse(startsWith(igraph::V(tree)$name, "node"), gsub(pattern = "node", replacement = "", igraph::V(tree)$name), igraph::V(tree)$name))
  }else{igraph::V(tree)$label <- NA}
  
  
  #3. Set node size
  if(is.null(node.size)){igraph::V(tree)$size <- 1}
  else if(node.size == "expansion"){
    # If node size is NULL (for germline), set node size to 1
    nodes <- lapply(nodes, function(x){
      if(is.null(x[["size"]])){x[["size"]] <- 1};return(x)
    })
    igraph::V(tree)$size <- as.numeric(lapply(nodes,function(x){x[["size"]]})[names(igraph::V(tree))])
  }
  # Scale the node size
  igraph::V(tree)$size <- igraph::V(tree)$size * node.size.scale
  
  #5. Set layout as tree and root on the germline
  layout <- igraph::layout_as_tree(tree, root = "germline")
  
  #4. Set edge length
  #igraph::E(tree)$length <- as.numeric(igraph::E(tree)$edge.length)
  
  # edge_df <- data.frame(length = as.numeric(igraph::E(tree)$edge.length), as_edgelist(tree))
  # layout_df <- data.frame(layout, node = V(tree)$name)
  # 
  # for (layer in )
  
  #6. Color the nodes
  # No node color based on features, use default colors and no pie charts
  if (is.null(node.color)){
    # Assign default node colors
    igraph::V(tree)$color <- ifelse(igraph::V(tree)$name == "germline", "orange",
                                    ifelse(startsWith(igraph::V(tree)$name, "node"), "lightblue", "grey"))
    
    #7. Plot tree
    igraph::plot.igraph(tree, layout = layout,
                        vertex.label = igraph::V(tree)$label,
                        edge.arrow.size = 0.1)
  }
  # Color on node features using pie charts
  else{
    #Create a list with the values to color per sequence per node
    color_values <- lapply(nodes,function(x){x[[node.color]]})[names(igraph::V(tree))]
    #Give the value "germline" to the germline node
    color_values$germline <- "germline"
    
    #Change NA values to "NA" for plotting purpose
    color_values <- lapply(color_values, function(x){
      index <- which(is.na(x))
      if (length(index) > 0){
        for (i in index){
          x[[i]] <- "NA"
        }
      }
      return(x)
    })
    
    #Get the unique values to color
    unique_values <- unique(unlist(color_values))
    
    #Specific color scheme for isotypes/VDJ_cgene
    if (node.color %in% c("isotype", "VDJ_cgene", "Isotype")){
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
    }
    #If there is no custom.color scheme provided, use default color scheme
    else if (is.null(custom.colors)){
      #Catch error of RColorBrewer when n < 3
      if (length(unique_values) > 2){color_list <- RColorBrewer::brewer.pal(length(unique_values), "Set1")}
      else if(length(unique_values) < 3){color_list <- RColorBrewer::brewer.pal(3, "Set1")}
    }
    #Use provided custom.colors
    else {color_list <- custom.colors}
    
    #Create list with the number of sequences per feature per node
    color_counts <- lapply(color_values, function(x){
      vector <- c()
      for (index in 1:length(unique_values)){
        count <- sum(stringr::str_count(x, pattern = unique_values[index]))
        vector[index] <- count
      }
      return(vector)
    })
    
    #7. Plot tree using pie charts
    igraph::plot.igraph(tree, layout = layout,
                        vertex.shape = "pie",
                        vertex.pie=color_counts,
                        vertex.pie.color=list(color_list),
                        vertex.label = igraph::V(tree)$label,
                        vertex.label.cex = 0.8,
                        edge.arrow.size = 0.1)
    #Add a legend
    graphics::legend(legend.position, legend = unique_values, fill = color_list, title = node.color,
                     cex = 0.7)      
        
  }
}
