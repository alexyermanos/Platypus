#' Plots lineage tree of clonotype from AntibodyForests object
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description ???
#' @param AntibodyForests_object AntibodyForests object - AntibodyForests object as obtained from the 'AntibodyForests()' function in Platypus.
#' @param sample string - denotes the sample that contains the clonotype.
#' @param clonotype string - denotes the clonotype from which the lineage tree should be plotted.
#' @param show.inner.nodes boolean - if TRUE, the tree with inner nodes is plotted (only present when the trees are created with the 'phylo.tree.nj', 'phylo.tree.mp', phylo.tree.ml', or 'phylo.tree.IgPhyML' construction algorithm). Defaults to FALSE.
#' @param color.by string - specifies the feature of the nodes that will be used for coloring the nodes. This sublist should be present in each sublist of each node in the 'nodes' objects within the AntibodyForests object. For each unique value for the selected feature, a unique color will be selected using the 'grDevices::rainbow()' function (unless a color gradient is created, see 'node.color.gradient' parameter). Defaults to NULL.
#' @param node.size string or integer or list of integers - specifies the size of the nodes. If set to 'expansion', the nodes will get a size that is equivalent to the number of cells that they represent. If set to an integer, all nodes will get this size. If set to a list of integers, in which each item is named to a node, the nodes will get these sizes. Defaults to 'expansion'.
#' @param node.size.factor integer - factor by which all node sizes are multiplied. Defaults to 1.
#' @param node.size.scale vector of integers - specifies the minimum and maximum node size in the plot, to which the number of cells will be scaled. Defaults to 10 and 30.
#' @param node.size.scale.range vector of integers - specifies the the range of the node size scale. Defaults to 1 and the first following multiple of 10 after the max node size.
#' @param node.color string or list of strings - specifies the color of nodes. If set to 'default', and the 'color.by' parameter is not specified, all the seqeuence-recovered nodes are colored lightblue. If set to 'default', and the 'color.by' parameter is set to a categorical value, the sequence-recovered nodes are colored  If set to a color (a color from the 'grDevices::color()' list or a valid HEX code), all the sequence-recovered nodes will get this color. If set to a list of colors, in which each item is named to a node, the nodes will get these colors. Defaults to 'default'.
#' @param node.color.gradient vector of strings - specifies the colors of the color gradient, if 'color.by' is set to a numerical feature. The minimum number of colors that need to be specified are 2. Defaults to 'c("#440154", "#481567", "#482677", "#453781", "#404788", "#39568C", "#33638D", "#2D708E", "#287D8E", "#238A8D", "#1F968B", "#20A387", "#29AF7F", "#3CBB75", "#55C667", "#73D055", "#95D840", "#B8DE29", "#DCE319", "#FDE725")'.
#' @param node.color.gradient.range vector of intgers - specifies the range of values over which the color gradient should be applied.
#' @param node.label string - specifies what should be plotted on the nodes. Options: 'name', 'size', and 'none'. Defaults to 'name'.
#' @param edge.label string - specifies what distance should be plotted on the edges. Options: 'dist' and 'none'. Defaults to 'dist'.
#' @param show.color.legend boolean - if TRUE, a legend is plotted to display the values of the specified node feature matched to the corresponding colors. Defaults to TRUE if the 'color.by' parameter is specified.
#' @param show.size.legend boolean - if TRUE, a legend is plotted to display the node sizes and the corresponding number of cells represented. Defaults to TRUE if the 'node.size' parameter is set to 'expansion'.
#' @param title string - specifies the title of the plot. Defaults to NULL.
#' @param color.legend.title string - specifies the title of the legend showing the color matching. Defaults to the name of the feature specified in the 'color.by' parameter.
#' @param size.legend.tile string - specifies the title of the legend showing the node sizes. Defaults to 'Expansion (# cells)'.
#' @return Plots lineage tree for the specified clonotype.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_plot(AntibodyForests_object,
#'                      sample = "S1",
#'                      clonotype = "clonotype1",
#'                      color.by = "Isotype",
#'                      main.title = "S1 clonotype1")
#'}

AntibodyForests_plot <- function(AntibodyForests_object,
                                 sample,
                                 clonotype,
                                 show.inner.nodes,
                                 color.by,
                                 node.size,
                                 node.size.factor,
                                 node.size.scale,
                                 node.size.scale.range,
                                 node.color,
                                 node.color.gradient,
                                 node.color.gradient.range,
                                 node.label,
                                 edge.label,
                                 show.color.legend,
                                 show.size.legend,
                                 main.title,
                                 color.legend.title,
                                 size.legend.title){
  
  predict_node_color_gradient_range <- function(values){
    
    # [function description]
    # Arguments:
    # - values:
    # Authors: Valentijn Tromp, Daphne van Ginneken
    
    #
    min_value <- min(values)
    max_value <- max(values)
    
    #
    nchar(strsplit(as.character(min_value), "\\.")[[1]][2])
    
    
  }
  
  predict_node_size_scale_range <- function(values){
    
    # [function description]
    # Arguments:
    # - values:
    # Authors: Valentijn Tromp, Daphne van Ginneken
    
    #
    min_value <- min(values)
    max_value <- max(values)
    
    #
    if(min_value < 10){floor_value <- 1}
    #
    else if (min_value < 100){floor_value <- floor(min_value/10)*10}
    #
    else{floor_value <- floor(min_value/10)*100}
    
    #
    if(max_value <= 10){ceiling_value <- 10}
    #
    else if(max_value <= 50){ceiling_value <- ceiling(max_value/10)*10}
    #
    else if(max_value <= 250){ceiling_value <- ceiling(max_value/50)*50}
    #
    else if(max_value <= 500){ceiling_value <- ceiling(max_value/100)*100}
    #
    else if(max_value <= 2500){ceiling_value <- ceiling(max_value/500)*500}
    #
    else{ceiling_value <- ceiling(max_value/1000)*1000}
    
    #
    return(c(floor_value, ceiling_value))
  }
  
  
  # 1. Retrieve objects from AntibodyForests object and perform input checks
  
  # If no AntibdoyForests object is provided, , a message is returned and execution is stopped
  if(missing(AntibodyForests_object)){stop("Please provide an AntibodyForests object that contains the lineage tree and its corresponding objects of the specified sample and clonotype!")}
  # If no sample and clonotype are specified, a message is returned and execution is stopped
  if(missing(sample) | missing(clonotype))stop("Please specify both a sample ID and clonotype ID, from the lineae of which the tree should be plotted!")
  
  # If the 'show.inner.nodes' parameter is not specified, it is set to FALSE
  if(missing(show.inner.nodes)){show.inner.nodes <- FALSE}
  
  # Retrieve igraph object from AntibodyForests object
  if(!show.inner.nodes){tree <- AntibodyForests_object[[sample]][[clonotype]][["igraph"]]}
  if(show.inner.nodes){tree <- AntibodyForests_object[[sample]][[clonotype]][["igraph.with.inner.nodes"]]}
  # If no tree could be found for the specified clonotype, a message is returned and execution is stopped
  if(is.null(tree) && !show.inner.nodes){stop(paste0("No tree could be found for ", clonotype, " of ", sample, "."))}
  if(is.null(tree) && show.inner.nodes){stop(paste0("No tree with inner nodes could be found for ", clonotype, " of ", sample, "."))}
  
  # If the feature specified to use for coloring the nodes could not be found for all nodes, a message is returned and execution is stopped
  if(!missing(color.by)){if(!(all(sapply(names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])[!names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) == "germline"], function(x) color.by %in% names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[x]]))))){stop("The feature specified with the 'color.by' parameter could not be found for all nodes.")}}
  # Retrieve the features (specified with the 'color.by' parameter) for each node from the AntibodyForests object and store the features in the 'node.feature.list', and replace all NA values with 'unknown'
  if(!missing(color.by)){node.feature.list <- as.list(sapply(names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])[!names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) == "germline"], function(x) AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[x]][[color.by]])); names(node.feature.list) <- names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])[!names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) == "germline"]; for(i in names(node.feature.list)){node.feature.list[[i]][is.na(node.feature.list[[i]])] <- "unknown"}}
  # If no feature is specified by the 'color.by' parameter, all the nodes will receive 'default' in the 'node.feature.list'
  if(missing(color.by)){node.feature.list <- as.list(sapply(names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])[!names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) == "germline"], function(x) "default")); names(node.feature.list) <- names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])[!names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) == "germline"]}
  
  # If the 'node.size' parameter is not specified, it is set to 'expansion'
  if(missing(node.size)){node.size <- "expansion"}
  # If the 'node.size' parameter is set to 'expansion', retrieve the size of each node from the AntibodyForests object and store the sizes in the 'node.size.list', and give all internal nodes and the germline node a size of 1
  if(is.character(node.size)){if(node.size == "expansion"){node.size.list <- as.list(sapply(names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])[!names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) == "germline"], function(x) AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[x]][["size"]])); names(node.size.list) <- names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])[!names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) == "germline"]; for(i in igraph::V(tree)$name){if(!i %in% names(node.size.list)){node.size.list[i] <- 1}}}}
  # If the 'node.size' parameter is set to a numerical value, all nodes will receive this size
  if(is.numeric(node.size)){node.size.list <- as.list(rep(node.size, length(igraph::V(tree)$name))); names(node.size.list) <- igraph::V(tree)$name}
  # If a list is provided with the 'node.size' parameter, this list is stored in the 'node.size.list' list
  if(is.list(node.size)){node.size.list <- node.size}
  # If the 'node.size.list' does not contain all the nodes, a message is returned and execution is stopped
  if(!all(igraph::V(tree)$name %in% names(node.size.list))){stop(paste(c("Not all nodes were found in the list specified by the 'node.size' parameter. Please include all node names:", paste(gtools::mixedsort(igraph::V(tree)$name), collapse = ", ")), collapse = " "))}
  # If the 'node.size.list' contains non-numerical or negative values, a message is returned and execution is stopped
  if(!(is.numeric(unlist(node.size.list)) && all(node.size.list > 0))){stop("Only positive numerical values are accepted by the 'node.size' parameter.")}
  
  # If the 'node.size.factor' parameter is not specified, it is set to 1
  if(missing(node.size.factor)){node.size.factor <- 1}
  # If the 'node.size.factor' parameter is set to a non-numerical or negative value, a message is returned and execution is stopped
  if(!(is.numeric(node.size.factor)) | !(if(is.numeric(node.size.factor)){node.size.factor > 0}else{FALSE})){stop("The 'node.size.factor' parameter only accepts positve numerical values.")}
  
  # If the 'node.size.scale.range' parameter is not specified, the minimum and maximum node size is set to the minimum and maximum value in the 'node.size.list' (multiplied by the 'node.size.factor')
  if(missing(node.size.scale.range)){node.size.scale.range <- predict_node_size_scale_range(unique(unlist(node.size.list)))}
  # If the 'node.size.scale.range' parameter contains non-numerical or negative values, a message is returned and execution is stopped
  if(!(is.numeric(node.size.scale.range) | !(if(is.numeric(node.size.scale.range)){all(node.size.scale.range > 0)}else{FALSE}) | length(node.size.scale.range) != 2)){stop("The 'node.size.scale.range' parameter only accepts a pair of positive numerical values.")}
  # If the 'node.size.scale.range' is set to a range that does not contain all values stored in the 'node.size.list' (multiplied by the 'node.size.factor'), a message is returned and execution is stopped
  if(!all(dplyr::between(c(min(unlist(node.size.list))*node.size.factor, max(unlist(node.size.list)))*node.size.factor, left = min(node.size.scale.range), right = max(node.size.scale.range)))){stop(paste0(c("The range specified with the 'node.size.scale.range' does not capture all the node sizes. The minimum range should be: ", min(unlist(node.size.list))*node.size.factor, " - ", max(unlist(node.size.list))*node.size.factor, ".")))}
  
  # If the 'node.size.scale' parameter is not specified, it is set to 'c(10, 20)'
  if(missing(node.size.scale)){node.size.scale <- c(10, 20)}
  # If the 'node.size.scale' parameter contains non-numerical or negative values, a message is returned and execution is stopped
  if(!(is.numeric(node.size.scale) | !(if(is.numeric(node.size.scale)){all(node.size.scale >= 0)}else{FALSE}) | length(node.size.scale) != 2)){stop("The 'node.size.scale' parameter only accepts a pair of positive numerical values.")}
  
  # Import the 'isotype_colors' list that specifies a unique color for each known isotype
  isotype_colors <- list(
    IgA   = "#8c96c6",  
    IgA1  = "#8856a7",  
    IgA2  = "#810f7c",  
    IgD   = "#f4a582",  
    IgE   = "#f7fcb9", 
    IgG   = "#238b45",
    IgG1  = "#41ab5d",
    IgG2  = "#74c476",
    IgG3  = "#a1d99b",
    IgG4  = "#c7e9c0",
    IgM   = "#084594",
    unknown = "#e0e0e0"
  )
  # Import the viridis color palette that can be used to create a color gradient
  viridis_palette <- c(
    "#440154",
    "#481567",
    "#482677",
    "#453781",
    "#404788",
    "#39568C",
    "#33638D",
    "#2D708E",
    "#287D8E",
    "#238A8D",
    "#1F968B",
    "#20A387",
    "#29AF7F",
    "#3CBB75",
    "#55C667",
    "#73D055",
    "#95D840",
    "#B8DE29",
    "#DCE319",
    "#FDE725"
  )
  
  # If the 'node.color' parameter is not specified, it is set to default
  if(missing(node.color)){node.color <- "default"}
  # If the 'node.color' parameter is set to 'default' and if only one unique value is present in the 'node.feature.list', the 'node.color' is set to 'lightblue'
  if(is.character(node.color)){if(node.color == "default" && length(unique(unlist(node.feature.list))) == 1){node.color <- "lightblue"}}
  # If the color specified by the 'node.color' parameter is not recognized, a message is returned and execution is stopped
  if(is.character(node.color)){if(node.color != "default" && !(node.color %in% grDevices::colors()) && !grepl(pattern = "^#([A-Fa-f0-9]{3}|[A-Fa-f0-9]{6})$", node.color)){stop("The color specified by the 'node.color' parameter is not recognized. Please choose a color from 'grDevices::colors()' or provide a hex code.")}}
  # If the 'node.color' parameter is set to a color, all nodes will be colors by this color
  if(is.character(node.color)){if(node.color != "default"){node.color.list <- as.list(rep(node.color, length(unique(unlist(node.feature.list))))); names(node.color.list) <- unique(unlist(node.feature.list))}}
  # If a list is provided with the 'node.color' parameter, this list is stored in the 'node.color.list' object
  if(is.list(node.color)){node.color.list <- node.color}
  # If no color gradient is specified with the 'node.color.gradient', and if not all nodes contain one unique numerical value for the selected feature in the 'node.feature.list', 'node.color.gradient' is set to 'none'
  if(missing(node.color.gradient)){if(!all(sapply(node.feature.list, function(x) length(unique(x)) == 1)) | !all(grepl(pattern = "^[0-9.-]+$", unique(unlist(node.feature.list))))){node.color.gradient <- "none"; node.color.gradient.range <- "none"; node.color.gradient.scale <- "none"}}
  # If no color gradient is specified with the 'node.color.gradient', and if all nodes contain only one unique numerical value for the selected feature in the 'node.feature.list', the 'viridis_palette' will be used to color the nodes
  if(missing(node.color.gradient)){if(all(sapply(node.feature.list, function(x) length(unique(x)) == 1)) && all(grepl(pattern = "^[0-9.-]+$", unique(unlist(node.feature.list))))){node.color.gradient <- viridis_palette}}
  # If not all nodes contain one unique numerical value for the selected feature in the 'node.feature.list', but the 'node.color.gradient' is not set to 'none', a message is returned and execution is stopped
  if(!(all(sapply(node.feature.list, function(x) length(unique(x)) == 1)) && all(grepl(pattern = "^[0-9.-]+$", unique(unlist(node.feature.list))))) && all(node.color.gradient != "none")){stop("The 'node.color.gradient' parameter can only be specified when only one unique numerical value is found per node for the feature selected in the 'color.by' parameter.")}
  # If the colors specified by the 'node.color.gradient' parameter are not all recognized, a message is returned and execution is stopped
  if(all(node.color.gradient != "none")){if(!all(node.color %in% grDevices::colors() | !grepl(pattern = "^#([A-Fa-f0-9]{3}|[A-Fa-f0-9]{6})$", node.color))){stop("Not all colors specified by the 'node.color.gradient' parameter are recognized. Please provide a vector with colors from 'grDevices::colors()' or with hex codes.")}}
  # If the 'node.color.gradient.range' is not specified, the min and max values in the 'node.feature.list' are used
  if(all(node.color.gradient != "none")){if(missing(node.color.gradient.range)){node.color.gradient.range <- c(min(as.numeric(unique(unlist(node.feature.list)))), max(as.numeric(unique(unlist(node.feature.list)))))}}
  # If the 'node.color.gradient.range' parameter does not contain a pair of numerical values, a message is returned and execution is stopped
  if(all(node.color.gradient != "none")){if(!(is.numeric(node.color.gradient.range) && length(node.color.gradient.range) == 2)){stop("The 'node.color.gradient.range' parameter only accepts a pair of numerical values.")}}
  # If the 'node.color.gradient.range' is set to a range that does not contain all values stored in the 'node.feature.list', a message is returned and execution is stopped
  if(all(node.color.gradient != "none")){if(!all(dplyr::between(c(min(as.numeric(unique(unlist(node.feature.list)))), max(as.numeric(unique(unlist(node.feature.list))))), left = min(node.color.gradient.range), right = max(node.color.gradient.range)))){stop(paste0(c("The range specified with the 'node.color.gradient.range' does not capture all the values of the feature selected in the 'color.by' parameter. The minimum range should be: ", min(as.numeric(unique(unlist(node.feature.list)))), " - ", max(as.numeric(unique(unlist(node.feature.list)))), ".")))}}
  # If 'node.color' is set to 'default', and all items in the 'node.feature.list' are known isotypes, the 'isotype_colors' list will be used to color the nodes 
  if(is.character(node.color)){if(node.color == "default" && all(node.color.gradient == "none") && all(unique(unlist(node.feature.list)) %in% names(isotype_colors))){node.color.list <- isotype_colors}}
  # If 'node.color' is set to 'default' and the 'node.color.gradient' is set to 'none', all unique values in the 'node.feature.list' will get a (random) color from the 'grDevices::rainbow()' function
  else if(is.character(node.color)){if(node.color == "default" && all(node.color.gradient == "none")){node.color.list <- as.list(grDevices::rainbow(length(unique(unlist(node.feature.list))))); names(node.color.list) <- sort(unique(unlist(node.feature.list)))}}
  # If 'node.color' is set to 'default' and the 'node.color.gradient' is not set to 'none', the specified color gradient will be used to assign a color to each numericla value that is present in the 'node.feature.list'
  if(class(node.color) == "character" && !missing(color.by)){if(node.color == "default" && all(node.color.gradient != "none")){
    # Retrieve the unique numerical values from the 'node.feature.list'
    numerical_values <- as.numeric(unique(unlist(node.feature.list)))
    # Calculate the largest number of digits after the decimal point, this 'factor' is used to determine how many colors need to be generated
    factor <- 10^max(sapply(numerical_values, function(x) nchar(strsplit(as.character(x), "\\.")[[1]][2])))
    # Calculate the number of colors that need to be generated by multiplying the absolute difference between the values in the 'node.color.gradient.range' vector by the 'factor' and substracting 1
    number_of_colors <- abs(max(node.color.gradient.range) - min(node.color.gradient.range))*factor-1
    # Create the color palette using the 'grDevices::colorRampPalette()' function
    color_palette <- grDevices::colorRampPalette(node.color.gradient)(number_of_colors)
    # To each unique numerical value, assign a color using the 'color_palette' vector
    node.color.list <- lapply(numerical_values, function(x) if(x == min(node.color.gradient.range)){return(color_palette[1])}else{color_palette[(x-min(node.color.gradient.range))*factor]}); names(node.color.list) <- unique(unlist(node.feature.list))
  }}
  # If the 'node.color.list' does not contain all the unique values in the 'node.feature.list', a message is returned and execution is stopped
  if(!all(unique(unlist(node.feature.list)) %in% names(node.color.list))){stop(paste(c("Not all values of the feature specified by the 'color.by' parameter are found in the list specified by the 'node.color' paramter. Please specify a color for the following values:", paste(gtools::mixedsort(unique(unlist(node.feature.list))), collapse = ", ")), collapse = " "))}
  # If the 'node.color.list' contains strings that are not recognized as a color, a message is returned and execution is stopped
  if(!all(sapply(unlist(node.color.list), function(x) x %in% grDevices::colors() | grepl(pattern = "^#([A-Fa-f0-9]{3}|[A-Fa-f0-9]{6})$", x)))){stop("Not all colors present in the list that is provided with the 'node.color' parameter are recognized.")}
  
  # If the 'node.label' parameter is not specified, it is set to 'name'
  if(missing(node.label)){node.label <- "name"}
  # If the 'node.label' parameter is not recognized, a message is returned and execution is stopped
  if(is.character(node.label)){if(!node.label %in% c("name", "size", "none")){stop("The 'node.label' parameter can only be set to 'name', 'size', 'none', or a list containing a size for all the nodes.")}}
  #
  if(is.list(node.label)){if(!all(igraph::V(tree)$name %in% names(node.label))){stop("")}}
  
  # If the 'show.edge.label' parameter is not specified, it is set to 'size', unless the 'show.inner.nodes' parameter is set to TRUE
  if(missing(edge.label) && show.inner.nodes){edge.label <- "none"}
  if(missing(edge.label) && !show.inner.nodes){edge.label <- "size"}
  
  # If the 'show.color.legend' parameter is not specified, it is set to FALSE, unless the 'color.by' parameter is specified
  if(missing(show.color.legend) && missing(color.by)){show.color.legend <- FALSE}
  if(missing(show.color.legend) && !missing(color.by)){show.color.legend <- TRUE}
  # If the 'show.size.legend' parameter is not specified, it is set to FALSE, unless the 'node.size' parameter is set to 'expansion'
  if(missing(show.size.legend) && node.size != "expansion"){show.size.legend <- FALSE}
  if(missing(show.size.legend) && node.size == "expansion"){show.size.legend <- TRUE}
  
  # If the 'title' parameter is not specified, it is set to an empty string
  if(missing(main.title)){main.title <- ""}
  # If the legend titles are not specified, the color legend title is set to the feature name and the size legend title is set to 'Expansion (# cells)'
  if(show.color.legend && missing(color.legend.title)){color.legend.title <- paste0(stringr::str_to_title(color.by))}
  if(show.size.legend && missing(size.legend.title)){size.legend.title <- "Expansion (# cells)"}
  
  
  # 2. Arrange the nodes in the igraph objects in a lineage tree format
  
  # Arrange the nodes by using 'germline' node as the root and by directing the tree downwards using the 'igraph::layout_as_tree()' function
  layout <- igraph::layout_as_tree(tree, root = "germline")
  
  
  # 3. Define the size of the nodes
  
  # Assign the node sizes from the 'node.size.list' to the igraph object
  igraph::V(tree)$size <- unlist(node.size.list[igraph::V(tree)$name])
  
  # Multiply the node sizes by the 'node.size.factor'
  igraph::V(tree)$size <- igraph::V(tree)$size*node.size.factor
  
  # Scale the node sizes using the 'node.size.scale' and 'node.size.scale.range' parameters (if the minimum size is not equal to the maximum size)
  if(min(igraph::V(tree)$size) != max(igraph::V(tree)$size)){igraph::V(tree)$size <- (igraph::V(tree)$size - min(node.size.scale.range)) / abs(diff(node.size.scale.range)) * abs(diff(node.size.scale)) + min(node.size.scale)}

  
  # 4. Define the color of the nodes
  
  # Iterate through the nodes
  for(i in 1:length(igraph::V(tree)$name)){
    # Retrieve the name of the current node
    node_name <- igraph::V(tree)$name[i]
    # If the name of the current node is 'germline', the node is plotted as a black circle
    if(node_name == "germline"){
      igraph::V(tree)$shape[i] = "circle"
      igraph::V(tree)$color[i] = "black"
      igraph::V(tree)$pie[i] = NA
      igraph::V(tree)$pie.color[i] = NA
    }
    # If the name of the current node starts with 'node' (which means it is a sequence-recovered node), the node is colored using the 'node.color' list
    else if(startsWith(node_name, prefix = "node")){
      # Retrieve the values for the feature (specified with the 'color.by' parameter) for the current node
      node_values <- node.feature.list[[node_name]]
      # If the current node contains only one unique value  for the specified feature, the node is plotted as a circle
      if(length(unique(node_values)) == 1){
        igraph::V(tree)$shape[i] = "circle"
        igraph::V(tree)$color[i] = node.color.list[[unique(node_values)]]
        igraph::V(tree)$pie[i] = NA
        igraph::V(tree)$pie.color[i] = NA
      }
      # If the current node contains multiple unique values for the specified feature, the node is plotted as a piechart 
      if(length(unique(node_values)) > 1){
        igraph::V(tree)$shape[i] = "pie"
        igraph::V(tree)$color[i] = NA
        igraph::V(tree)$pie[i] = list(as.vector(unlist(lapply(names(node.color.list), function(x) sum(node_values == x)))))
        igraph::V(tree)$pie.color[i] = list(as.vector(unlist(node.color.list)))
      }
    }
    # If the name of the current node is not 'germline' and does not start with 'node', it is a unrecovered sequence node, that will be colored grey
    else{
      igraph::V(tree)$shape[i] = "circle"
      igraph::V(tree)$color[i] = "white"
      igraph::V(tree)$pie[i] = NA
      igraph::V(tree)$pie.color[i] = NA
    }
  }
  
  
  # 5. Define the node and edge labels
  
  # If 'node.label' is set to 'name', the germline node becomes 'G', the sequence-recovered nodes only keep their node number, and the remaining nodes will loose their label
  if(is.character(node.label)){if(node.label == "name"){
    igraph::V(tree)$label <- ifelse(igraph::V(tree)$name == "germline", 
                                    yes = "G", 
                                    no = ifelse(startsWith(igraph::V(tree)$name, "node"), 
                                                yes = gsub(pattern = "node", replacement = "", igraph::V(tree)$name), 
                                                no = ""))
  }}
  
  # If 'node.label' is set to 'size', the original number of cells represent by the nodes are displayed as node labels
  if(is.character(node.label)){if(node.label == "size"){
    igraph::V(tree)$label <- sapply(igraph::V(tree)$name, function(x) if(x == "germline"){return("G")}else{return(as.character(node.size.list[[as.character(x)]]))})
  }}
  
  # If a list is provided as 'node.label' object, the labels of this list are assigned to the nodes
  if(is.list(node.label)){igraph::V(tree)$label <- node.label[igraph::V(tree)$name]}
  
  # If 'node.label' is set to 'none', the node labels are set to NA before plotting the tree
  if(is.character(node.label)){if(node.label == "none"){igraph::V(tree)$label <- NA}}
  
  # For all black nodes, change the color of the label to 'white'
  igraph::V(tree)$label.color[igraph::V(tree)$color == "black"] <- "white"
  
  # Resize the size of the node labels according to the size of the node itself
  igraph::V(tree)$label.cex <- igraph::V(tree)$size/15
  
  
  # 6. Plot lineage tree
  
  # If a title is specified, add 2 lines of margin on top of the plot, and if the 'show.color.legend' or 'show.size.legend' is set to TRUE, add 6 lines of margin on the right side of the plot
  if(main.title == "" && !show.color.legend && !show.size.legend){par(mar = c(0, 0, 0, 0), xpd = TRUE)}
  if(main.title != "" && !show.color.legend && !show.size.legend){par(mar = c(0, 0, 2, 0), xpd = TRUE)}
  if(main.title == "" && (show.color.legend | show.size.legend)){par(mar = c(0, 0, 0, 6), xpd = TRUE)}
  if(main.title != "" && (show.color.legend | show.size.legend)){par(mar = c(0, 0, 2, 6), xpd = TRUE)}
  
  # Plot the tree using the 'igraph::plot.igraph()' function 
  igraph::plot.igraph(tree,
                      layout = layout,
                      vertex.label.dist = 0,
                      edge.arrow.size = 0.25)
  
  
  # 7. Add legend(s)
  
  # If 'show.color.legend' and 'show.size.legend' are both set to TRUE, the color legend is positioned above, while the size legend is positioned below
  if(show.color.legend && show.size.legend){
    y_color_legend <- 0.9*par("usr")[4]
    y_size_legend <- 0 -0.2*par("usr")[4]
  }
  # If only the 'show.color.legend' is set to TRUE, the color legend is positioned above
  else if(show.color.legend && !show.size.legend){
    y_color_legend <- 0.9*par("usr")[4]
    y_size_legend <- NA
  }
  # If only the 'size.color.legend' is set to TRUE, the size legend is positioned above
  else if(!show.color.legend && show.size.legend){
    y_color_legend <- NA
    y_size_legend <- 0.9*par("usr")[4]
  }
  
  # If 'show.color.legend' is set to TRUE, add a legend to the plot to show the color matching
  if(show.color.legend){
    # If 'node.color.gradient' is set to 'none', a legend is created showing the matching of the colors with the unique values found for the feature selected in the 'color.by' parameter 
    if(all(node.color.gradient == "none")){
      #
      legend_values <- c(gtools::mixedsort(unique(unlist(node.feature.list))))
      #
      legend_colors <- c(unlist(node.color.list[gtools::mixedsort(unique(unlist(node.feature.list)))]))
      #
      color_legend <- graphics::legend(legend = legend_values,       #
                                       pt.bg = legend_colors,        #
                                       pch = 21,                     #
                                       pt.cex = 1.5,                 #
                                       title = color.legend.title,   #
                                       title.adj = 0,                #
                                       bty = "n",                    #
                                       x = par("usr")[2],            #
                                       y = y_color_legend,           #
                                       cex = 0.75)                   #
    }
    # If node colors are specified with a gradient, create a legend using a color gradient
    if(all(node.color.gradient != "none")){
      #
      gradient_values <- c(min(node.color.gradient.range), rep(NA, 1000), mean(c(min(node.color.gradient.range), mean(node.color.gradient.range))), rep(NA, 1000), mean(node.color.gradient.range), rep(NA, 1000), mean(c(mean(node.color.gradient.range), max(node.color.gradient.range))), rep(NA, 1000), max(node.color.gradient.range))
      #
      gradient_values[!is.na(gradient_values)] <- round(gradient_values[!is.na(gradient_values)], 3)
      #
      gradient_values[!is.na(gradient_values)] <- format(gradient_values[!is.na(gradient_values)], 3)
      #
      gradient_colors <- grDevices::colorRampPalette(node.color.gradient)(4005)
      #
      color_legend <- graphics::legend(legend = rev(gradient_values),         #
                                       fill = rev(gradient_colors),           #
                                       border = NA,                           #
                                       y.intersp = c(1, rep(0.0025, 4003)),   #
                                       title = color.legend.title,            #
                                       title.adj = 0,                         #
                                       bty = "n",                             #
                                       x = par("usr")[2],                     #
                                       y = y_color_legend,                    #
                                       cex = 0.75)                            #
    }
  }
  
  # If 'show.size.legend' is set to TRUE, add a legend to the plot to show the node sizes
  if(show.size.legend){
    #
    legend_values <- c(min(node.size.scale.range), mean(node.size.scale.range)-0.5, max(node.size.scale.range))
    #
    legend_values <- round(legend_values, 0)
    #
    if(max(node.size.scale) > 20){legend_values <- paste(paste(rep(" ", max(node.size.scale)-20)/5, collapse = ""), legend_values)}
    #
    legend_sizes <- c(min(node.size.scale), mean(node.size.scale), max(node.size.scale))/200
    #
    size_legend <- graphics::legend(legend = legend_values,                                      #
                                    pt.cex = 0,                                                  #
                                    y.intersp = c(1.25, rep(max(node.size.scale)/10+0.25, 2)),   #
                                    title = size.legend.title,                                   #
                                    title.adj = 0,                                               #
                                    bty = "n",                                                   #
                                    x = par("usr")[2],                                           #
                                    y = y_size_legend,                                           #
                                    cex = 0.75)                                                  #
    #
    x_positions <- (size_legend$text$x + size_legend$rect$left) / 2
    y_positions <- size_legend$text$y
    #
    graphics::symbols(x = x_positions, y = y_positions, circles = legend_sizes, inches = FALSE, add = TRUE, bg = 'black')
  }
  
  
  # 8. Add title to the plot
  if(main.title != ""){mtext(main.title, font = 2)}
}