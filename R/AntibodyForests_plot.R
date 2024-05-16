#' Plots lineage tree of clonotype from AntibodyForests object
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description This function retrieves the igraph object from the provided AntibodyForests object for the specified clone within the specified sample and plots the lineage tree using the specified plotting parameters.
#' @param AntibodyForests_object AntibodyForests object - AntibodyForests object as obtained from the 'AntibodyForests()' function in Platypus.
#' @param sample string - denotes the sample that contains the clonotype.
#' @param clonotype string - denotes the clonotype from which the lineage tree should be plotted.
#' @param show.inner.nodes boolean - if TRUE, the tree with inner nodes is plotted (only present when the trees are created with the 'phylo.tree.nj', 'phylo.tree.mp', phylo.tree.ml', or 'phylo.tree.IgPhyML' construction algorithm). Defaults to FALSE.
#' @param color.by string - specifies the feature of the nodes that will be used for coloring the nodes. This sublist should be present in each sublist of each node in the 'nodes' objects within the AntibodyForests object. For each unique value for the selected feature, a unique color will be selected using the 'grDevices::rainbow()' function (unless a color gradient is created, see 'node.color.gradient' parameter). Defaults to NULL.
#' @param label.by string - specifies what should be plotted on the nodes. Options: 'name', 'size', a feature that is stored in the 'nodes' list, and 'none'. Defaults to 'name'.
#' @param edge.label string - specifies what distance between the nodes is shown as labels of the edges. Options: 'original' (distance that is stored in the igraph object), 'none' (no edge labels are shown), 'lv' (Levensthein distance), 'dl' (Damerau-Levenshtein distance), 'osa' (Optimal String Alignment distance), and 'hamming' (Hamming distance). Defaults to 'lv'. 
#' @param node.size string or integer or list of integers - specifies the size of the nodes. If set to 'expansion', the nodes will get a size that is equivalent to the number of cells that they represent. If set to an integer, all nodes will get this size. If set to a list of integers, in which each item is named to a node, the nodes will get these sizes. Defaults to 'expansion'.
#' @param node.size.factor integer - factor by which all node sizes are multiplied. Defaults to 1.
#' @param node.size.scale vector of integers - specifies the minimum and maximum node size in the plot, to which the number of cells will be scaled. Defaults to 10 and 30.
#' @param node.size.range vector of integers - specifies the the range of the node size scale. Defaults to 1 and the first following multiple of 10 after the max node size.
#' @param node.color string or list of strings - specifies the color of nodes. If set to 'default', and the 'color.by' parameter is not specified, all the seqeuence-recovered nodes are colored lightblue. If set to 'default', and the 'color.by' parameter is set to a categorical value, the sequence-recovered nodes are colored  If set to a color (a color from the 'grDevices::color()' list or a valid HEX code), all the sequence-recovered nodes will get this color. If set to a list of colors, in which each item is named to a node, the nodes will get these colors. Defaults to 'default'.
#' @param node.color.gradient vector of strings - specifies the colors of the color gradient, if 'color.by' is set to a numerical feature. The minimum number of colors that need to be specified are 2. Defaults to 'c("#440154", "#481567", "#482677", "#453781", "#404788", "#39568C", "#33638D", "#2D708E", "#287D8E", "#238A8D", "#1F968B", "#20A387", "#29AF7F", "#3CBB75", "#55C667", "#73D055", "#95D840", "#B8DE29", "#DCE319", "#FDE725")'.
#' @param show.color.legend boolean - if TRUE, a legend is plotted to display the values of the specified node feature matched to the corresponding colors. Defaults to TRUE if the 'color.by' parameter is specified.
#' @param show.size.legend boolean - if TRUE, a legend is plotted to display the node sizes and the corresponding number of cells represented. Defaults to TRUE if the 'node.size' parameter is set to 'expansion'.
#' @param main.title string - specifies the main title of the plot. Defaults to NULL.
#' @param color.legend.title string - specifies the title of the legend showing the color matching. Defaults to the (capitalized) name of the feature specified in the 'color.by' parameter (converted by the 'stringr::str_to_title()' function).
#' @param size.legend.tile string - specifies the title of the legend showing the node sizes. Defaults to 'Expansion (# cells)'.
#' @return Plots lineage tree for the specified clonotype.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_plot(AntibodyForests_object,
#'                      sample = "S1",
#'                      clonotype = "clonotype1",
#'                      color.by = "isotype",
#'                      main.title = "Lineage tree of S1 clonotype1")
#'}

AntibodyForests_plot <- function(AntibodyForests_object,
                                 sample,
                                 clonotype,
                                 show.inner.nodes,
                                 color.by,
                                 label.by,
                                 node.size,
                                 node.size.factor,
                                 node.size.scale,
                                 node.size.range,
                                 node.color,
                                 node.color.gradient,
                                 edge.label,
                                 show.color.legend,
                                 show.size.legend,
                                 main.title,
                                 color.legend.title,
                                 size.legend.title){
  
  predict_range_labels <- function(values){
    
    # Determines the most optimal set of numbers to describe a range, whereby the set is returned as a vector of two floats in scientific notation
    # Arguments:
    # - values: list of numerical values (integers and/or floats)
    # Authors: Valentijn Tromp, Daphne van Ginneken
    
    # Convert values to numeric
    values <- as.numeric(values)
    
    # Find the minimum and maximum value in the list
    min_value <- min(values)
    max_value <- max(values)
    
    # Determine the order of magnitude for the minimum and maximum values
    min_value_order_of_magnitude <- floor(log10(abs(min_value)))
    max_value_order_of_magnitude <- floor(log10(abs(max_value)))
    
    # Convert the minimum and maximum value to the scientific notation
    min_value <- paste(as.character(round(min_value/10^min_value_order_of_magnitude, digits = 1)), "e", sprintf(fmt = "%+03d", min_value_order_of_magnitude), sep = "")
    max_value <- paste(as.character(round(max_value/10^max_value_order_of_magnitude, digits = 1)), "e", sprintf(fmt = "%+03d", max_value_order_of_magnitude), sep = "")
    
    # Format minimum and maximum values by making sure there is one digit after the decimal point
    min_value <- gsub(pattern = ".*e", replacement = paste0(sprintf(fmt = "%.1f", as.numeric(stringr::str_extract(pattern = ".*(?=e)", min_value))), "e"), min_value)
    max_value <- gsub(pattern = ".*e", replacement = paste0(sprintf(fmt = "%.1f", as.numeric(stringr::str_extract(pattern = ".*(?=e)", max_value))), "e"), max_value)
    
    # If the maximum value is rounded to a multiple of 10, correct the scientific notation
    if(gsub(pattern = "\\..*", replacement = "", max_value) == 10){
      max_value <- gsub(pattern = "10\\.", replacement = "1.", max_value)
      max_value <- gsub(pattern = paste0("\\", sprintf(fmt = "%+03d", max_value_order_of_magnitude)), replacement = sprintf(fmt = "%+03d", max_value_order_of_magnitude+1), max_value)
    }
    if(gsub(pattern = "\\..*", replacement = "", max_value) == -10){
      max_value <- gsub(pattern = "-10\\.", replacement = "-1.", max_value)
      max_value <- gsub(pattern = paste0("\\", sprintf(fmt = "%+03d", max_value_order_of_magnitude)), replacement = sprintf(fmt = "%+03d", max_value_order_of_magnitude-1), max_value)
    }
    
    # If the minimum value is rounded to a multiple of 10, correct the scientific notation
    if(gsub(pattern = "\\..*", replacement = "", min_value) == 10){
      min_value <- gsub(pattern = "10\\.", replacement = "1.", min_value)
      min_value <- gsub(pattern = paste0("\\", sprintf(fmt = "%+03d", min_value_order_of_magnitude)), replacement = sprintf(fmt = "%+03d", min_value_order_of_magnitude+1), min_value)
    }
    if(gsub(pattern = "\\..*", replacement = "", min_value) == -10){
      min_value <- gsub(pattern = "-10\\.", replacement = "-1.", min_value)
      min_value <- gsub(pattern = paste0("\\", sprintf(fmt = "%+03d", min_value_order_of_magnitude)), replacement = sprintf(fmt = "%+03d", min_value_order_of_magnitude-1), min_value)
    }
    
    # Extract the numbers before and after the decimal point (necessary to determine the appropriate range)
    min_value_number_before_point <- as.numeric(gsub(pattern = "\\..*", replacement = "", min_value))
    min_value_number_behind_point <- as.numeric(gsub(pattern = "(.*\\.)|(e.*)", replacement = "", min_value))
    max_value_number_before_point <- as.numeric(gsub(pattern = "\\..*", replacement = "", max_value))
    max_value_number_behind_point <- as.numeric(gsub(pattern = "(.*\\.)|(e.*)", replacement = "", max_value))
    
    # If the minimum value is positive...
    if(as.numeric(min_value) > 0){
      # and if the number behind the decimal point is between 0 and 5, replace this number by 5
      if(min_value_number_behind_point > 0 && min_value_number_behind_point < 5){
        min_value <- base::gsub(pattern = "\\.\\de", replacement = ".5e", min_value)
      }
      # and if the number behind the decimal point is between 5 and 10, replace this number by 0 (and add +1 to the number before the decimal point)
      if(min_value_number_behind_point > 5 && min_value_number_behind_point < 10){
        min_value <- base::gsub(pattern = ".*\\.*", replacement = paste0(as.character(as.numeric(min_value_number_before_point)+1), "."), min_value)
        min_value <- base::gsub(pattern = "\\.\\de", replacement = ".0e", min_value)
      }
      # and if 9 became 10 before the decimal point, correct the scientific notation
      if(gsub(pattern = "\\..*", replacement = "", max_value) == 10){
        min_value <- gsub(pattern = "10\\.", replacement = "1.", min_value)
        min_value <- gsub(pattern = paste0("\\", sprintf(fmt = "%+03d", min_value_order_of_magnitude)), replacement = sprintf(fmt = "%+03d", min_value_order_of_magnitude+1), min_value)
      }
    }
    
    # If the minimum value is negative...
    if(as.numeric(min_value) < 0){
      # and if the number behind the decimal point is between 0 and 5, replace this number by 0
      if(min_value_number_behind_point > 0 && min_value_number_behind_point < 5){
        min_value <- base::gsub(pattern = "\\.\\de", replacement = ".0e", min_value)
      }
      # and if the number behind the decimal point is between 5 and 10, replace this number by 5 
      if(min_value_number_behind_point > 5 && min_value_number_behind_point < 10){
        min_value <- base::gsub(pattern = "\\.\\de", replacement = ".5e", min_value)
      }
    }
    
    # If the maximum value is positive...
    if(as.numeric(max_value) > 0){
      # and if the number behind the decimal point is between 0 and 5, replace this number by 0
      if(max_value_number_behind_point > 0 && max_value_number_behind_point < 5){
        max_value <- base::gsub(pattern = "\\.\\de", replacement = ".0e", max_value)
      }
      # and if the number behind the decimal point is between 0 and 10, replace this number by 5 
      if(max_value_number_behind_point > 5 && max_value_number_behind_point < 10){
        max_value <- base::gsub(pattern = "\\.\\de", replacement = ".5e", max_value)
      }
    }
    
    # If the maximum value is negative...
    if(as.numeric(max_value) < 0){
      # and if the number behind the decimal point is between 0 and 5, replace this number by 5
      if(max_value_number_behind_point > 0 && max_value_number_behind_point < 5){
        max_value <- base::gsub(pattern = "\\.\\de.", replacement = ".5e", max_value)
      }
      # and if the number behind the decimal point is between 5 and 10, replace this number by 0 (and substract -1 frpm the number before the decimal point)
      if(max_value_number_behind_point > 5 && max_value_number_behind_point < 10){
        max_value <- base::gsub(pattern = ".*\\.", replacement = paste0(as.character(as.numeric(max_value_number_before_point)-1), "."), max_value)
        max_value <- base::gsub(pattern = "\\.\\de", replacement = ".0e", max_value)
      }
      # and if -9 became -10 before the decimal point, correct the scientific notation
      if(gsub(pattern = "\\..*", replacement = "", max_value) == -10){
        max_value <- gsub(pattern = "-10\\.", replacement = "-1.", max_value)
        max_value <- gsub(pattern = paste0("\\", sprintf(fmt = "%+03d", max_value_order_of_magnitude)), replacement = sprintf(fmt = "%+03d", max_value_order_of_magnitude-1), max_value)
      }
    }
    
    # Return the range
    return(c(as.numeric(min_value), as.numeric(max_value)))
  }
  
  predict_node_size_scale_range <- function(values){
    
    # Determines the most optimal range for a list of node sizes, whereby the range is returned as a vector of two integers 
    # Arguments:
    # - values: list of integers
    # Authors: Valentijn Tromp, Daphne van Ginneken
    
    # Convert values to numeric
    values <- as.numeric(values)
    
    # Find the minimum and maximum value in the list
    min_value <- min(values)
    max_value <- max(values)
    
    # If the minimum value is below 10, set the minimum value of the range to 1
    if(min_value < 10){min_value <- 1}
    # Else, set the minimum value of the range to the nearest lower multiple of the logarithm of the minimum node size
    else{min_value <- floor(min_value/10^floor(log10(min_value)))*10^floor(log10(min_value))}
    
    # If the maximum value is equal to or below 10, set the maximum value of the range to 10
    if(max_value <= 10){min_value <- 10}
    # Else, set the maximum value of the range to the nearest higher multiple of the logarithm of the maximum node size
    else{max_value <- ceiling(max_value/10^floor(log10(max_value)))*10^floor(log10(max_value))}
    
    # If the difference between the maximum value of the range and the original maximum value is bigger than 0.5 of the logarithm of the maximum node size, subtract this from the maximum value
    if((max_value - max(values)) > 0.5*10^floor(log10(max_value))){max_value <- max_value - 0.5*10^floor(log10(max_value))}
    
    # Return the range
    return(c(min_value, max_value))
  }
  
  
  # 1. Retrieve objects from AntibodyForests object and perform input checks
  
  # If no AntibdoyForests object is provided, a message is returned and execution is stopped
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
  
  # If the 'label.by' parameter is not specified, it is set to 'name'
  if(missing(label.by)){label.by <- "name"}
  # If the 'label.by' parameter is not recognized, a message is returned and execution is stopped
  if(!missing(label.by)){if(!label.by %in% c("name", "size", "none")){if(!(all(sapply(names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])[!names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) == "germline"], function(x) label.by %in% names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[x]]))))){stop("The feature specified with the 'label.by' parameter could not be found for all nodes.")}}}
  # Retrieve the features (specified with the 'label.by' parameter) for each node from the AntibodyForests object and store the features in the 'node.feature.list', and replace all NA values with 'unknown'
  if(!missing(color.by)){if(!label.by %in% c("name", "size", "none")){node.label.list <- as.list(sapply(names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])[!names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) == "germline"], function(x) unique(AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[x]][[label.by]]))); names(node.label.list) <- names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])[!names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) == "germline"]; for(i in names(node.label.list)){node.label.list[[i]][is.na(node.label.list[[i]])] <- "unknown"}}}
  
  # If the 'edge.label' parameter is not specified, and inner nodes are plotted, it is set to 'original'
  if(missing(edge.label) && show.inner.nodes){edge.label <- "original"}
  # If the 'edge.label' parameter is not specified, and no inner nodes are plotted, it is set to 'lv' (Levensthein distance)
  if(missing(edge.label) && !show.inner.nodes){edge.label <- "lv"}
  # If the 'edge.label' parameter is not recognized, a message is returned and execution is stopped
  if(!edge.label %in% c("original", "none", "lv", "dl", "osa", "hamming")){stop("The 'edge.label' parameter is not recognized. Please choose from the following options: 'original', 'none', 'lv', 'dl', 'osa', 'hamming'.")}
  # If the 'edge.label' is set to string distance metric, while 'show.inner.nodes' is set to TRUE, a message is returned and execution is stopped 
  if(show.inner.nodes && !edge.label %in% c("original", "none")){stop("When non-recovered internal nodes are present in the tree, the 'edge.label' parameter can only be set to 'original' (to show the distance that is stored in the igraph object), or 'none'.")}
  
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
  
  # If the 'node.size.range' parameter is not specified, the scale is determined using the 'predict_node_size_scale()' function
  if(missing(node.size.range)){node.size.range <- predict_node_size_scale_range(unique(unlist(node.size.list)))}
  # If the 'node.size.range' parameter contains non-numerical or negative values, a message is returned and execution is stopped
  if(!(is.numeric(node.size.range) | !(if(is.numeric(node.size.range)){all(node.size.range > 0)}else{FALSE}) | length(node.size.range) != 2)){stop("The 'node.size.range' parameter only accepts a pair of positive numerical values.")}
  # If the 'node.size.range' is set to a range that does not contain all values stored in the 'node.size.list' (multiplied by the 'node.size.factor'), a message is returned and execution is stopped
  if(!all(dplyr::between(c(min(unlist(node.size.list))*node.size.factor, max(unlist(node.size.list)))*node.size.factor, left = min(node.size.range), right = max(node.size.range)))){stop(paste0(c("The range specified with the 'node.size.range' does not capture all the node sizes. The minimum range should be: ", min(unlist(node.size.list))*node.size.factor, " - ", max(unlist(node.size.list))*node.size.factor, ".")))}
  
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
  # Import the viridis color palette from the 'scales::viridis_pal()' function that can be used to create a color gradient with the 'scales::pal_seq_gradient()' function
  viridis_palette <- scales::viridis_pal()(1000)
  
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
  
  # If no color gradient is specified with the 'node.color.gradient' parameter, and if not all nodes contain one unique numerical value for the selected feature in the 'node.feature.list', 'node.color.gradient' is set to 'none'
  if(missing(node.color.gradient)){if(!all(sapply(node.feature.list, function(x) length(unique(x)) == 1)) | !all(grepl(pattern = "^[0-9.-]+$", unique(unlist(node.feature.list))))){node.color.gradient <- "none"}}
  # If no color gradient is specified with the 'node.color.gradient' parameter, and if all nodes contain only one unique numerical value for the selected feature in the 'node.feature.list', the 'viridis_palette' will be used to color the nodes
  if(missing(node.color.gradient)){if(all(sapply(node.feature.list, function(x) length(unique(x)) == 1)) && all(grepl(pattern = "^[0-9.-]+$", unique(unlist(node.feature.list))))){node.color.gradient <- viridis_palette}}
  # If both a color and color gradient are specified, a message is returned and execution is stopped
  if(is.character(node.color)){if(node.color != "default" && !all(node.color.gradient == "none")){stop("Both the 'node.color' and the 'node.color.gradient' are specified. Please specify one.")}}
  # If not all nodes contain one unique numerical value for the selected feature in the 'node.feature.list', but the 'node.color.gradient' is not set to 'none', a message is returned and execution is stopped
  if(!(all(sapply(node.feature.list, function(x) length(unique(x)) == 1)) && all(grepl(pattern = "^[0-9.-]+$", unique(unlist(node.feature.list))))) && all(node.color.gradient != "none")){stop("The 'node.color.gradient' parameter can only be specified when only one unique numerical value is found per node for the feature selected in the 'color.by' parameter.")}
  # If the colors specified by the 'node.color.gradient' parameter are not all recognized, a message is returned and execution is stopped
  if(all(node.color.gradient != "none")){if(!all(node.color %in% grDevices::colors() | !grepl(pattern = "^#([A-Fa-f0-9]{3}|[A-Fa-f0-9]{6})$", node.color))){stop("Not all colors specified by the 'node.color.gradient' parameter are recognized. Please provide a vector with colors from 'grDevices::colors()' or with hex codes.")}}
  
  # If 'node.color' is set to 'default', and all items in the 'node.feature.list' are known isotypes, the 'isotype_colors' list will be used to color the nodes 
  if(is.character(node.color)){if(node.color == "default" && all(node.color.gradient == "none") && all(unique(unlist(node.feature.list)) %in% names(isotype_colors))){node.color.list <- isotype_colors}}
  # If 'node.color' is set to 'default' and the 'node.color.gradient' is set to 'none', all unique values in the 'node.feature.list' will get a (random) color from the 'grDevices::rainbow()' function
  else if(is.character(node.color)){if(node.color == "default" && all(node.color.gradient == "none")){node.color.list <- as.list(grDevices::rainbow(length(unique(unlist(node.feature.list))))); names(node.color.list) <- sort(unique(unlist(node.feature.list)))}}
  # If 'node.color' is set to 'default' and the 'node.color.gradient' is not set to 'none', the specified color gradient will be used to assign a color to each numericla value that is present in the 'node.feature.list'
  if(class(node.color) == "character" && !missing(color.by)){if(node.color == "default" && all(node.color.gradient != "none")){
    # For each node, retrieve the unique numerical values from the 'node.feature.list'
    numerical_values <- as.numeric(unique(unlist(node.feature.list)))
    # Assign a color to each value/node using the 'scale::cscale()' function from a color gradient that is created using the 'scales::pal_seq_gradient()' function with the colors specified in the 'node.color.gradient' vector 
    node.color.list <- as.list(scales::cscale(numerical_values, palette = scales::pal_gradient_n(node.color.gradient)))
    # Set the names of the colors in the 'node.color.lst' to the numerical values to which they are matched
    names(node.color.list) <- as.character(numerical_values)
    }}
  # If the 'node.color.list' does not contain all the unique values in the 'node.feature.list', a message is returned and execution is stopped
  if(!all(unique(unlist(node.feature.list)) %in% names(node.color.list))){stop(paste(c("Not all values of the feature specified by the 'color.by' parameter are found in the list specified by the 'node.color' parameter. Please specify a color for the following values:", paste(gtools::mixedsort(unique(unlist(node.feature.list))), collapse = ", ")), collapse = " "))}
  # If the 'node.color.list' contains strings that are not recognized as a color, a message is returned and execution is stopped
  if(!all(sapply(unlist(node.color.list), function(x) x %in% grDevices::colors() | grepl(pattern = "^#([A-Fa-f0-9]{3}|[A-Fa-f0-9]{6})$", x)))){stop("Not all colors present in the list that is provided with the 'node.color' parameter are recognized.")}
  
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
  igraph::V(tree)$size <- sapply(igraph::V(tree)$name, function(x){
    if(x == "germline"){return(1)}
    else if(startsWith(x, "node")){return(node.size.list[[x]])}
    else{return(0.5)}
  })
  
  # Multiply the node sizes by the 'node.size.factor'
  igraph::V(tree)$size <- igraph::V(tree)$size*node.size.factor
  
  # Scale the node sizes using the 'node.size.scale' and 'node.size.range' parameters (if the minimum size is not equal to the maximum size)
  if(min(igraph::V(tree)$size) != max(igraph::V(tree)$size)){
    igraph::V(tree)$size <- sapply(igraph::V(tree)$name, function(x){
      if(x == "germline"){return(min(node.size.scale))}
      else if(startsWith(x, "node")){return((node.size.list[[x]] - min(node.size.range)) / abs(diff(node.size.range)) * abs(diff(node.size.scale)) + min(node.size.scale))}
      else{return(0.5*min(node.size.scale))}
    })
  }
  
  
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
  
  # If 'label.by' is set to 'name', the germline node becomes 'G', the sequence-recovered nodes only keep their node number, and the remaining nodes will loose their label
  if(label.by == "name"){
    igraph::V(tree)$label <- ifelse(igraph::V(tree)$name == "germline", 
                                    yes = "G", 
                                    no = ifelse(startsWith(igraph::V(tree)$name, "node"), 
                                                yes = gsub(pattern = "node", replacement = "", igraph::V(tree)$name), 
                                                no = ""))
  }
  
  # If 'label.by' is set to 'size', the original number of cells represent by the nodes are displayed as node labels
  if(label.by == "size"){
    igraph::V(tree)$label <- sapply(igraph::V(tree)$name, function(x){
      if(x == "germline"){return("G")}
      else if(startsWith(x, "node")){return(as.character(AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[x]][["size"]]))}
      else{return("")}
    })
  }
  
  # If a list is provided as 'node.label' object, the labels of this list are assigned to the nodes
  if(!label.by %in% c("name", "size", "none")){
    igraph::V(tree)$label <- sapply(igraph::V(tree)$name, function(x){
      if(x == "germline"){return("G")}
      else if(startsWith(x, "node")){return(as.character(node.label.list[[x]]))}
      else{return("")}
    })
  }
  
  # If the 'label.by' parameter is set to 'none', the node labels are set to NA before plotting the tree
  if(label.by == "none"){igraph::V(tree)$label <- NA}
  
  # For all black nodes, change the color of the label to 'white'
  igraph::V(tree)$label.color[igraph::V(tree)$color == "black"] <- "white"
  
  # Resize the size of the node labels according to the size of the node itself
  igraph::V(tree)$label.cex <- igraph::V(tree)$size/10
  
  # If the 'edge.label' is set to 'default', display the edge lengths from the graphs as node labels in the graph
  if(edge.label == "original"){igraph::E(tree)$label <- as.character(round(as.numeric(igraph::E(tree)$edge.length), 1))}
  # If the 'edge.label' parameter is set to a string distance metric, calculate the distance between the nodes with the 'stringdist::stringdist()' function
  if(edge.label %in% c("osa", "lv", "dl", "hamming")){
    # Retrieve all the items that are stored for each node in the 'nodes' list in the AntibodyForest object
    node_items <- unique(unlist(sapply(AntibodyForests_object[[sample]][[clonotype]][["nodes"]], function(x) names(x))))
    # Select the names of the sequences
    seq_names <- unique(unlist(sapply(node_items, function(x) if(all(sapply(AntibodyForests_object[[sample]][[clonotype]][["nodes"]], function(y) is.character(y[[x]]) && length(y[[x]]) == 1))){return(x)})))
    # Assign a label to each edge 
    igraph::E(tree)$label <- sapply(igraph::as_ids(igraph::E(tree)), function(x){
      # Retrieve the nodes that are connected by this edge
      nodes <- strsplit(x, split = "\\|")[[1]]
      # Calculate the distance between these nodes by summing the distance for all sequences
      dist <- sum(sapply(seq_names, function(y){
        # Retrieve this sequence for both nodes 
        seq1 <- AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[(nodes[1])]][[y]]
        seq2 <- AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[(nodes[2])]][[y]]
        # Return the distance between nodes used the specified distance metric
        return(stringdist::stringdist(seq1, seq2, method = edge.label))
      }))
    })
  }
  # If the 'edge.label' parameter is set to "none", the edge labels are set to NA before plotting the tree
  if(edge.label == "none"){igraph::E(tree)$label <- NA}
  
  
  # 6. Plot lineage tree
  
  # Add 2 lines of margin to all sides of the pland, and if the 'show.color.legend' or 'show.size.legend' is set to TRUE, add 6 lines of margin on the right side of the plot
  if(show.color.legend | show.size.legend){par(mar = c(2, 2, 2, 6), xpd = TRUE)}
  else{par(mar = c(2, 2, 2, 2), xpd = TRUE)}
  
  # Plot the tree using the 'igraph::plot.igraph()' function 
  igraph::plot.igraph(tree,                     # Plot the tree and its attributes (if present)
                      layout = layout,          # Place the vertices on the plot using the layout of a lineage tree
                      vertex.label.dist = 0,    # Center the node labels on the nodes
                      edge.arrow.size = 0.25)   # Set the size of the arrows
  
  
  # 7. Add legend(s)
  
  # If 'show.color.legend' is set to TRUE, add a legend to the plot to show the color matching
  if(show.color.legend){
    # If 'node.color.gradient' is set to 'none', a legend is created showing the matching of the colors with the unique values found for the feature selected in the 'color.by' parameter 
    if(all(node.color.gradient == "none")){
      # Retrieve the values of the feature and sort them
      legend_values <- c(gtools::mixedsort(unique(unlist(node.feature.list))))
      # Assign the corresponding colors to the values 
      legend_colors <- c(unlist(node.color.list[gtools::mixedsort(unique(unlist(node.feature.list)))]))
      # Create the color legend (without a gradient)
      color_legend <- graphics::legend(legend = legend_values,       # Assign legend values to the legend points
                                       pt.bg = legend_colors,        # Specify the colors for the legend points
                                       pch = 21,                     # Use filled circles as legend points
                                       pt.cex = 1.5,                 # Set the size of the legend circles
                                       title = color.legend.title,   # Specify the title of the color legend
                                       title.cex = 0.8,              # Set the size of the legend title
                                       title.adj = 0.25,             # Align the title
                                       bty = "n",                    # Do not draw a box around the legend
                                       x = par("usr")[2],            # Set x-position of the legend (right side of the graph)
                                       y = par("usr")[4]*0.8,        # Set y-position of the legend (top of the graph)
                                       cex = 0.6)                    # Set the size of the legend
    }
    # If node colors are specified with a gradient, create a legend using a color gradient
    if(all(node.color.gradient != "none")){
      # Create the colors of the gradient with the 'scales::pal_seq_gradient()' function
      gradient_colors <- grDevices::colorRampPalette(node.color.gradient)(1000)
      # Calculate the minimum and maximum label of the color gradient with the 'predict_range_labels()' function
      gradient_text_labels <- predict_range_labels(numerical_values)
      # Calculate the three values that lie in between the minimum and maximum label
      gradient_text_labels <- c(min(gradient_text_labels), mean(c(min(gradient_text_labels), mean(gradient_text_labels))), mean(gradient_text_labels), mean(c(mean(gradient_text_labels), max(gradient_text_labels))), max(gradient_text_labels))
      # If all absolute values fall within the range of 0.1 to 10, scientific notation is not utilized
      if(abs(min(gradient_text_labels)) > 0.1 && abs(max(gradient_text_labels)) < 100){
        gradient_text_labels <- format(gradient_text_labels, scientific = FALSE, digits = 3)
      } 
      # Else, values are not expressed in scientific notation
      else{
        gradient_text_labels <- format(gradient_text_labels, scientific = TRUE, digits = 3)
      }
      # Calculate a vector in which the five values are matched to the colors in the 'gradient_colors' vector by placing the values at the right position in the vector 
      gradient_labels <- rep(" ", 1000)
      # Iterate through the five labels that will be displayed next to the color gradient
      for(i in gradient_text_labels){
        # If the label equals to the minimum value of the range, the position of the label is set to 1 (on the bottom)
        if(as.numeric(i) == min(numerical_values)){
          position <- 1
          # At this position in the 'gradient_labels' vector (modified by adding +30 in order to be at the right position next to the color gradient), add the label
          gradient_labels[position+30] <- i; 
          # At the three surrounding positions in the 'gradient_colors' vector, change the color to white
          gradient_colors[position:(position+3)] <- rep("#FFFFFF", 3)
          }
        # If the label does not equal to the minimum value, the position of the label needs to be calculated
        else{
          position <- round((as.numeric(i)-min(numerical_values))/abs(min(numerical_values) - max(numerical_values))*1000)
          # At this position in the 'gradient_labels' vector (modified by adding +30 in order to be at the right position next to the color gradient), add the label
          gradient_labels[position+30] <- i
          # At the five surrounding positions in the 'gradient_colors' vector, change the color to white
          gradient_colors[(position-1):(position+1)] <- rep("#FFFFFF", 3)
          }
        }
      # Create the color legend (with a gradient)
      color_legend <- graphics::legend(legend = rev(gradient_labels),             # Assign gradient values to the color gradient in reverse order 
                                       fill = rev(gradient_colors),               # Fill the legend points (boxes) with the gradien colors in reverse order
                                       border = NA,                               # Do not draw a border around the legend points
                                       y.intersp = c(0.5, rep(0.01, 999)),        # Set the vertical spacing between legend points (which creates the gradient)
                                       title = paste(color.legend.title, "\n"),   # Specify the title of the color legend
                                       title.cex = 0.8,                           # Set the size of the legend title
                                       title.adj = 0.25,                          # Align the title 
                                       bty = "n",                                 # Do not draw a box around the legend
                                       x = par("usr")[2],                         # Set x-position of the legend (right side of the graph)
                                       y = par("usr")[4]*0.8,                     # Set y-position of the legend (top of the graph)
                                       cex = 0.6)                                 # Set the size of the legend
    }
  }
  
  # If a color legend is plotted, position the size legend underneath it  
  if(show.color.legend){
    y_size_legend <-  par("usr")[4]*0.8 - 1.05*color_legend$rect$h
  } 
  # Else, position the size legend in the upper right corner
  else{
    y_size_legend <-  par("usr")[4]*0.8
  }
  
  # If 'show.size.legend' is set to TRUE, add a legend to the plot to show the node sizes
  if(show.size.legend){
    # Determine the three node sizes that will be displayed in the size legend, which include the min, mean, and max value in the 'node.size.range'
    legend_values <- c(min(node.size.range), mean(node.size.range)-0.5, max(node.size.range))
    # Round the legend values
    legend_values <- round(legend_values, 0)
    # Assign sizes for the legend (dividing the sizes by 200 assures that the node sizes in the legend correspond to the node sizes in the graph, and this factor is retrieved from the source code of the 'igraph::plot.igraph()' function)
    legend_sizes <- c(min(node.size.scale), mean(node.size.scale), max(node.size.scale))/200
    # Create the size legend
    size_legend <- graphics::legend(legend = legend_values,                                                               # Assign legend values to the nodes
                                    pt.cex = 0.6,                                                                         # Mask legend points by setting size to zero
                                    y.intersp = c(1, 0.6*max(node.size.scale)/10+0.1, 0.9*max(node.size.scale)/10+0.1),   # Set vertical spacing between the text of the legend
                                    title = size.legend.title,                                                            # Specify the title of the size legend
                                    title.cex = 0.8,                                                                      # Set the size of the legend title
                                    title.adj = 0.25,                                                                     # Align the title 
                                    bty = "n",                                                                            # Do not draw a box around the legend
                                    x = par("usr")[2],                                                                    # Set x-position of the legend (right side of the graph)
                                    y = y_size_legend,                                                                    # Set y-position of the legend (top of the graph or below size legend)
                                    cex = 0.8)                                                                            # Set the size of the legend
    # Determine the position of the circles in the plot
    x_positions <- (size_legend$text$x + size_legend$rect$left) / 2 
    y_positions <- size_legend$text$y
    # Add the circles to the plot
    graphics::symbols(x = x_positions, y = y_positions, circles = legend_sizes, inches = FALSE, add = TRUE, bg = 'black')
  }
  
  
  # 8. Add title to the plot
  if(main.title != ""){mtext(main.title, font = 2)}
}