#' Plots lineage tree of clonotype from AntibodyForests object
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description This function retrieves the igraph object from the provided AntibodyForests object for the specified clone within the specified sample and plots the lineage tree using the specified plotting parameters.
#' @param AntibodyForests_object AntibodyForests object - AntibodyForests object as obtained from the 'AntibodyForests()' function in Platypus.
#' @param sample string - denotes the sample that contains the clonotype.
#' @param clonotype string - denotes the clonotype from which the lineage tree should be plotted.
#' @param show.inner.nodes boolean - if TRUE, the tree with inner nodes is plotted (only present when the trees are created with the 'phylo.tree.nj', 'phylo.tree.mp', phylo.tree.ml', or 'phylo.tree.IgPhyML' construction algorithm). Defaults to FALSE.
#' @param x.scaling
#' @param y.scaling
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
#' @param sub.title string - specifies the sub title of the plot that is to be plotted below the main title. Defaults to NULL.
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
                                 x.scaling,
                                 y.scaling,
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
                                 sub.title,
                                 color.legend.title,
                                 size.legend.title){
  
  
  calculate_optimal_x_scaling <- function(tree){
    
    #
    node_names <- igraph::V(tree)$name
    
    #
    node_sizes <- igraph::V(tree)$size
    names(node_sizes) <- node_names
    
    #
    node_coordinates <- as.data.frame(igraph::norm_coords(igraph:::i.postprocess.layout(igraph::layout_as_tree(tree, root = "germline")), xmin = -1, xmax = 1, ymin = -1, ymax = 1))
    rownames(node_coordinates) <- node_names
    colnames(node_coordinates) <- c("x", "y")
    
    #
    min_distances <- c()
    min_distances_names <- c()
    
    # Loop over all the y coordinates
    for(i in unique(node_coordinates$y)){
      # Create an empty vector to store the horizontal distances between the nodes at this height
      horizontal_distances <- c()
      # Select the coordinates of the nodes at the current height
      current_nodes <- node_coordinates[node_coordinates$y == i, ]
      # If there are multiple nodes present at this height in the lineage tree...
      if(nrow(current_nodes) > 1){
        # Loop over these nodes
        for(j in rownames(current_nodes)){
          # For each possible pair of nodes at this height...
          for(k in rownames(current_nodes)[rownames(current_nodes) != j]){
            # Calculate the distance between the nodes by calculating the absolute difference between the x coordinates and subsequently substracting the node sizes
            dist <- abs(node_coordinates[j, "x"] - node_coordinates[k, "x"])  - node_sizes[[j]]/200 - node_sizes[[k]]/200
            names(dist) <- paste(j, k, sep = "-")
            # Append this distance to the 'horizontal_distances' vector
            horizontal_distances <- c(horizontal_distances, dist) 
          }
        }
        # Select the minimum distance from the 'horizontal_distances' vector 
        min_distance <- min(horizontal_distances)
        min_distance_name <- names(horizontal_distances[horizontal_distances == min_distance])[1]
        # Add this distance to the 'min_distances' vector
        min_distances <- c(min_distances, min_distance)
        min_distances_names <- c(min_distances_names, min_distance_name)
        names(min_distances) <- min_distances_names
      }
    }
    
    #
    if(length(min_distances) != 0){
      #
      nodes_with_biggest_overlap <- names(min_distances)[min_distances == min(min_distances)][1]
      #
      node_a <- strsplit(nodes_with_biggest_overlap, split = "-")[[1]][1]
      node_b <- strsplit(nodes_with_biggest_overlap, split = "-")[[1]][2]
      #
      current_distance <- abs(node_coordinates[node_a, "x"] - node_coordinates[node_b, "x"])
      #
      theoretical_distance <- node_sizes[[node_a]]/200 + node_sizes[[node_b]]/200 + 0.2*min(node_sizes)/200
      #
      x_scaling_factor <- theoretical_distance/current_distance
      #
      if(x_scaling_factor < 1){x_scaling_factor <- 1}
      #
      return(x_scaling_factor)
    }
    
    #
    if(length(min_distances) == 0){
      return(1)
    }
  }
  
  
  plot_igraph_object <- function(x,
                                 xlim, 
                                 ylim,
                                 legend,
                                 title,
                                 ...){
    
    # Save input object as 'graph' and ensure that the input object is an igraph object
    graph <- x
    igraph:::ensure_igraph(graph)
    # Count the number of vertices in the graph
    vc <- igraph::vcount(graph)
    # Parse through the plot parameters 
    params <- igraph:::i.parse.plot.params(graph, list(...))
    # Retrieve the node size (and divide by 200) and shape
    vertex.size <- 1/200 * params("vertex", "size")
    shape <- igraph:::igraph.check.shapes(params("vertex", "shape"))
    # Retrieve the color palette to use for vertex color (if provided)
    palette <- params("plot", "palette")
    if(!is.null(palette)){old_palette <- palette(palette); on.exit(palette(old_palette), add = TRUE)}
    # Retrieve the node label properties
    label.family <- params("vertex", "label.family")
    label.font <- params("vertex", "label.font")
    label.cex <- params("vertex", "label.cex")
    label.degree <- params("vertex", "label.degree")
    label.color <- params("vertex", "label.color")
    label.dist <- params("vertex", "label.dist")
    labels <- params("vertex", "label")
    # Retrieve the edge properties
    edge.color <- params("edge", "color")
    edge.width <- params("edge", "width")
    edge.lty <- params("edge", "lty")
    arrow.size <- params("edge", "arrow.size")[1]
    arrow.width <- params("edge", "arrow.width")[1]
    edge.labels <- params("edge", "label")
    loop.angle <- params("edge", "loop.angle")
    # Retrieve and process the arrows (if present) and their mode/direction
    arrow.mode <- igraph:::i.get.arrow.mode(graph, params("edge", "arrow.mode"))
    # Retrieve the edge curvature (if provided)
    curved <- params("edge", "curved")
    if(is.function(curved)){curved <- curved(graph)}
    # Retrieve the edge label properties
    edge.label.font <- params("edge", "label.font")
    edge.label.family <- params("edge", "label.family")
    edge.label.cex <- params("edge", "label.cex")
    edge.label.color <- params("edge", "label.color")
    elab.x <- params("edge", "label.x")
    elab.y <- params("edge", "label.y")
    # Retrieve and process the tree layout
    layout <- igraph:::i.postprocess.layout(params("plot", "layout"))
    # Retrieve other plot parameters (margin, scaling, aspect ratio, and frame)
    margin <- params("plot", "margin")
    margin <- rep(margin, length.out = 4)
    rescale <- params("plot", "rescale")
    asp <- params("plot", "asp")
    frame.plot <- params("plot", "frame.plot")
    # Retrieve titles
    main <- params("plot", "main")
    sub <- params("plot", "sub")
    xlab <- params("plot", "xlab")
    ylab <- params("plot", "ylab")
    # Retrieve the maximum vertex size
    maxv <- max(vertex.size)
    # If 'rescale' is set to TRUE, rescale the layout to the specified 'xlim' and 'ylim'
    if(rescale){
      layout <- igraph::norm_coords(layout, xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2])
    }
    # If a legend is to be plotted on the right side of the plot, add +1.5 to the upper limit of the horizontal axis
    if(legend){xlim[2] <- xlim[2]+1.5}
    # If a title is to be plotted above the plot, add +0.5 to the upper limit of the vertical axis
    if(title){ylim[2] <- ylim[2]+0.5}
    # Set up the plot 
    graphics::plot(x = 0, y = 0, type = "n",   # Center the plot and produce no points or lines (yet)
                   xlim = xlim, ylim = ylim,   # Set the x and y limits
                   axes = FALSE,               # Draw no axes on the plot
                   frame.plot = FALSE,         # Draw no box around the plot
                   asp = asp,                  # Set the y/x aspect ratio
                   main = main, sub = sub,     # Specicy the main and sub title of the plot
                   xlab = xlab, ylab = ylab)   # Specify the axis titles
    # Retrieve a list of the edges in the graph
    el <- igraph::as_edgelist(graph, names = FALSE)
    # Retrieve the edge labels 
    edge.labels <- edge.labels
    # Create matrix to store edge coordinatores
    edge.coords <- matrix(0, nrow = nrow(el), ncol = 4)
    # Prepare edge coordinates
    edge.coords[, 1] <- layout[, 1][el[, 1]]   # x0 (x top)
    edge.coords[, 2] <- layout[, 2][el[, 1]]   # y0 (y top)
    edge.coords[, 3] <- layout[, 1][el[, 2]]   # x1 (x down)
    edge.coords[, 4] <- layout[, 2][el[, 2]]   # y1 (y down)
    # If all the nodes have the same shape, clip both ends of the edges to the nodes using the function stored in the 'igraph:::.igraph.shapes' object
    if(length(unique(shape)) == 1){
      ec <- igraph:::.igraph.shapes[[shape[1]]]$clip(edge.coords, el, params = params, end = "both")
    }
    # If not, clip the ends of the edges separately using the different functions stored in the 'igraph:::.igraph.shapes' object
    else{
      shape <- rep(shape, length.out = igraph::vcount(graph))
      ec <- edge.coords
      ec[, 1:2] <- t(sapply(seq(length.out = nrow(el)), function(x){
        igraph:::.igraph.shapes[[shape[el[x, 1]]]]$clip(edge.coords[x, , drop = FALSE], el[x, , drop = FALSE], params = params, end = "from")
      }))
      ec[, 3:4] <- t(sapply(seq(length.out = nrow(el)), function(x){
        igraph:::.igraph.shapes[[shape[el[x, 2]]]]$clip(edge.coords[x, , drop = FALSE], el[x, , drop = FALSE], params = params, end = "to")
      }))
    }
    # Store the edge coordinates in separate objects ('x0', 'y0', 'x1', and 'y1')
    x0 <- ec[, 1]
    y0 <- ec[, 2]
    x1 <- ec[, 3]
    y1 <- ec[, 4]
    # Plot edges with the appropriate properties
    if(length(edge.color) > 1){edge.color <- edge.color}
    if(length(edge.width) > 1){edge.width <- edge.width}
    if(length(edge.lty) > 1){edge.lty <- edge.lty}
    if(length(arrow.mode) > 1){arrow.mode <- arrow.mode}
    if(length(arrow.size) > 1){arrow.size <- arrow.size}
    if(length(curved) > 1){curved <- curved}
    # If the same arrow is to be used for all the edges, plot the edges as this type of arrow with the appropriate properties using the 'igraph:::igraph.Arrows()' function
    if(length(unique(arrow.mode)) == 1){
      lc <- igraph:::igraph.Arrows(x0, y0, x1, y1, 
                                   h.col = edge.color, sh.col = edge.color, 
                                   h.lwd = 1, sh.lwd = edge.width, 
                                   open = FALSE, 
                                   code = arrow.mode[1], 
                                   sh.lty = edge.lty, h.lty = 1, 
                                   size = arrow.size, width = arrow.width,
                                   curved = curved)
      lc.x <- lc$lab.x
      lc.y <- lc$lab.y
    }
    # If different type of arrows are to be used, plot the arrows separately
    else{
      curved <- rep(curved, length.out = ecount(graph))
      lc.x <- lc.y <- numeric(length(curved))
      for(code in 0:3){
        valid <- arrow.mode == code
        if(!any(valid)){next}
        ec <- edge.color
        if(length(ec) > 1){ec <- ec[valid]}
        ew <- edge.width
        if(length(ew) > 1){ew <- ew[valid]}
        el <- edge.lty
        if(length(el) > 1){el <- el[valid]}
        lc <- igraph.Arrows(x0[valid], y0[valid], x1[valid], y1[valid],
                            code = code, 
                            sh.col = ec, h.col = ec, 
                            sh.lwd = ew, h.lwd = 1, 
                            h.lty = 1, sh.lty = el, 
                            open = FALSE,
                            size = arrow.size, width = arrow.width, 
                            curved = curved[valid])
        lc.x[valid] <- lc$lab.x
        lc.y[valid] <- lc$lab.y
      }
      # Adjust edge label positions (if specified)
      if(!is.null(elab.x)){lc.x <- ifelse(is.na(elab.x), lc.x, elab.x)}
      if(!is.null(elab.y)) {lc.y <- ifelse(is.na(elab.y), lc.y, elab.y)}
      # Plot the edge labels
      text(x = lc.x, y = lc.y, labels = edge.labels, col = edge.label.color, family = edge.label.family, font = edge.label.font, cex = edge.label.cex)
    }
    # If all the nodes have the same shape, plot the nodes with the appropriate properties using the function stored in the 'igraph:::.igraph.shapes' object
    if(length(unique(shape)) == 1){igraph:::.igraph.shapes[[shape[1]]]$plot(layout, params = params)}
    # If not, plot the nodes differently using the functions stored in the 'igraph:::.igraph.shapes' object
    else{sapply(seq(length.out = igraph::vcount(graph)), function(x){igraph:::.igraph.shapes[[shape[x]]]$plot(layout[x, , drop = FALSE], v = x, params = params)})}
    # Restore the graphical parameters after function execution
    old_xpd <- par(xpd = TRUE)
    on.exit(par(old_xpd), add = TRUE)
    # Define node label positions
    x <- layout[, 1]
    y <- layout[, 2] 
    # If all the node labels have the same font family, plot all the node labels using this font family
    if(length(label.family) == 1){text(x, y, labels = labels, col = label.color, family = label.family, font = label.font, cex = label.cex)}
    # If the node labels have different font families, apply these fonts separately
    else{
      if1 <- function(vect, idx){
        if(length(vect) == 1){vect}
        else{vect[idx]}
      }
      sapply(seq_len(igraph::vcount(graph)), function(v){text(x[v], y[v], labels = if1(labels, v), col = if1(label.color, v), family = if1(label.family, v), font = if1(label.font, v), cex = if1(label.cex, v))
      })}
  }
  
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
  
  # If the 'x.scaling' and 'y.sclaing' parameters are not defined, they are set to 'default'
  if(missing(x.scaling)){x.scaling <- "default"}
  if(missing(y.scaling)){y.scaling <- "default"}
  
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
  if(is.numeric(node.size)){node.size.list <- as.list(rep(node.size, length(names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])))); names(node.size.list) <- names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])}
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
  
  #
  if(missing(node.size.range) && length(unique(node.size.list)) == 1){node.size.range <- c(unlist(unique(node.size.list)), unlist(unique(node.size.list)))}
  # If the 'node.size.range' parameter is not specified and the 'node.size.list' contains different sizes, the scale is determined using the 'predict_node_size_scale()' function
  if(missing(node.size.range) && length(unique(node.size.list)) > 1){node.size.range <- predict_node_size_scale_range(unique(unlist(node.size.list)))}
  # If the 'node.size.range' parameter contains non-numerical or negative values, a message is returned and execution is stopped
  if(!(is.numeric(node.size.range) | !(if(is.numeric(node.size.range)){all(node.size.range > 0)}else{FALSE}) | length(node.size.range) != 2)){stop("The 'node.size.range' parameter only accepts a pair of positive numerical values.")}
  # If the 'node.size.range' is set to a range that does not contain all values stored in the 'node.size.list' (multiplied by the 'node.size.factor'), a message is returned and execution is stopped
  if(!all(dplyr::between(c(min(unlist(node.size.list))*node.size.factor, max(unlist(node.size.list)))*node.size.factor, left = min(node.size.range), right = max(node.size.range)))){stop(paste0(c("The range specified with the 'node.size.range' does not capture all the node sizes. The minimum range should be: ", min(unlist(node.size.list))*node.size.factor, " - ", max(unlist(node.size.list))*node.size.factor, ".")))}
  
  #
  if(missing(node.size.scale) && length(unique(node.size.list)) == 1){node.size.scale <- c(unlist(unique(node.size.list)), unlist(unique(node.size.list)))}
  # If the 'node.size.scale' parameter is not specified and the 'node.size.list' contains different sizes, it is set to 'c(10, 20)'
  if(missing(node.size.scale) && length(unique(node.size.list)) > 1){node.size.scale <- c(10, 20)}
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
  
  # If the 'main.title' or 'sub.title' parameter are not specified, they are set to an empty string
  if(missing(main.title)){main.title <- ""}
  if(missing(sub.title)){sub.title <- ""}
  
  # If the legend titles are not specified, the color legend title is set to the feature name and the size legend title is set to 'Expansion (# cells)'
  if(show.color.legend && missing(color.legend.title)){color.legend.title <- paste0(stringr::str_to_title(color.by))}
  if(show.size.legend && missing(size.legend.title)){size.legend.title <- "Expansion (# cells)"}
  
  
  # 2. Arrange the nodes in the igraph objects in a lineage tree format
  
  # Arrange the nodes by using 'germline' node as the root and by directing the tree downwards using the 'igraph::layout_as_tree()' function
  layout <- igraph::layout_as_tree(tree, root = "germline")
  
  
  # 3. Define the size of the nodes
  
  # Assign the node sizes from the 'node.size.list' to the igraph object
  #igraph::V(tree)$size <- sapply(igraph::V(tree)$name, function(x){
  #  if(x == "germline"){return(1)}
  #  else if(startsWith(x, "node")){return(node.size.list[[x]])}
  #  else{return(0.5)}
  #})
  igraph::V(tree)$size <- sapply(igraph::V(tree)$name, function(x) node.size.list[[x]])
  
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
  igraph::V(tree)$label.cex <- 0.75
  
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
  
  
  # 6. Define horizontal and vertical scaling
  
  #
  if(all(x.scaling == "default")){
    x.scaling <- c(-calculate_optimal_x_scaling(tree), calculate_optimal_x_scaling(tree))
  }
  
  #
  if(all(y.scaling == "default")){
    y.scaling <- c((-length(unique(layout[, 2]))-1)/8, (length(unique(layout[, 2]))-1)/8)
    if(max(y.scaling) < 1){y.scaling <- c(-1, 1)}
    if(max(y.scaling) < 0.5*max(x.scaling)){y.scaling <- 0.5*x.scaling}
  }
  
  
  # 7. Plot lineage tree
  
  # Plot the tree using the 'igraph::plot.igraph()' function 
  plot_igraph_object(tree,                                              # Plot the tree and its attributes (if present)
                     layout = layout,                                   # Place the vertices on the plot using the layout of a lineage tree
                     xlim = x.scaling,                                  # Set the limits for the horizontal axis, thereby scaling the horizontal node distance
                     ylim = y.scaling,                                  # Set the limits for the vertical axis, thereby scaling the vertical node distance
                     legend = (show.color.legend | show.size.legend),   # Specify whether empty space should be plotted on the right side 
                     title = if(main.title != ""){TRUE}else{FALSE},     # Specify whether empty space should be plotted on top 
                     edge.arrow.size = 0.1)                             # Set the size of the arrows
  
  
  # 8. Add legend(s)
  
  # If 'show.size.legend' is set to TRUE, add a legend to the plot to show the node sizes
  if(show.size.legend){
    # Determine the three node sizes that will be displayed in the size legend, which include the min, mean, and max value in the 'node.size.range'
    legend_values <- c(min(node.size.range), mean(node.size.range)-0.5, max(node.size.range))
    # Round the legend values
    legend_values <- round(legend_values, 0)
    #
    legend_values <- paste0(legend_values, paste(rep(" ", 50), collapse = ""))
    # Assign sizes for the legend (dividing the sizes by 200 assures that the node sizes in the legend correspond to the node sizes in the graph, and this factor is retrieved from the source code of the 'igraph::plot.igraph()' function, see line 10)
    legend_node_sizes <- c(min(node.size.scale), mean(node.size.scale), max(node.size.scale))/200
    #
    x_positions <- rep(x.scaling[2]+0.2+max(legend_node_sizes), 3)
    #
    y_positions <- c(-0.25-legend_node_sizes[1], -0.25-legend_node_sizes[1]*2-0.05-legend_node_sizes[2], -0.25-legend_node_sizes[1]*2-0.05-legend_node_sizes[2]*2-0.05-legend_node_sizes[3])
    #
    graphics::text(x = x.scaling[2]+0.2,
                   y = -0.1,
                   label = size.legend.title,
                   adj = 0)
    #
    graphics::symbols(x = x_positions,
                      y = y_positions,
                      circles = legend_node_sizes, 
                      inches = FALSE,
                      bg = "black",
                      add = TRUE)
    #
    graphics::text(x = x_positions + max(legend_node_sizes) + 0.1,
                   y = y_positions,
                   labels = legend_values,
                   adj = 0,
                   cex = 0.75)
  }
  
  # If 'show.color.legend' is set to TRUE, add a legend to the plot to show the color matching
  if(show.color.legend){
    # If 'node.color.gradient' is set to 'none', a legend is created showing the matching of the colors with the unique values found for the feature selected in the 'color.by' parameter 
    if(all(node.color.gradient == "none")){
      # Retrieve the values of the feature and sort them
      legend_values <- c(gtools::mixedsort(unique(unlist(node.feature.list))))
      # Assign the corresponding colors to the values 
      legend_colors <- c(unlist(node.color.list[gtools::mixedsort(unique(unlist(node.feature.list)))]))
      #
      legend_circle_size <- 1.0 / length(legend_values) / 2 - 0.025
      #
      if(legend_circle_size > 0.1){legend_circle_size <- 0.1}
      #
      if(show.size.legend){x_positions <- rep((x.scaling[2]+0.2+max(legend_node_sizes)), length(legend_values))}
      if(!show.size.legend){x_positions <- rep((x.scaling[2]+0.2+legend_circle_size), length(legend_values))}
      #
      y_positions <- sapply(1:length(legend_values), function(x) (y.scaling[2]-0.1) - (x*2-1)*(0.1+0.025))
      #
      graphics::text(x = x.scaling[2]+0.2,
                     y = y.scaling[2]-0.1,
                     label = color.legend.title,
                     adj = 0)
      #
      graphics::symbols(x = x_positions,
                        y = y_positions,
                        circles = rep(legend_circle_size, length(legend_values)), 
                        inches = FALSE,
                        bg = legend_colors,
                        add = TRUE)
      #
      graphics::text(x = x_positions + legend_circle_size + 0.05,
                     y = y_positions,
                     labels = legend_values,
                     adj = 0,
                     cex = 0.75)
      
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
      #
      if(show.size.legend){x_positions <- rep((x.scaling[2]+0.2+max(legend_node_sizes)), 1000)}
      if(!show.size.legend){x_positions <- rep((x.scaling[2]+0.2+0.1), 1000)}
      #
      y_positions <- rev(seq(0.25, y.scaling[2]-0.25, (abs(y.scaling[2])-0.5)/999))
      #
      graphics::text(x = x.scaling[2]+0.2,
                     y = y.scaling[2]-0.1,
                     label = color.legend.title,
                     adj = 0)
      #
      graphics::symbols(x = x_positions,
                        y = y_positions,
                        rectangles = matrix(rep(c(0.2, (y.scaling[2]-0.25)/1000), 1000), ncol = 2, nrow = 1000, byrow = TRUE),
                        lty = 0, 
                        inches = FALSE,
                        bg = rev(gradient_colors),
                        add = TRUE)
      #
      graphics::symbols(x = c(rep(unique(x_positions)+0.075, 5), rep(unique(x_positions)-0.075, 5)),
                        y = rep(sapply(gradient_text_labels, function(x) return(0.25+round((as.numeric(x)-min(numerical_values))/abs(min(numerical_values) - max(numerical_values))*1000)*(abs(y.scaling[2])-0.5)/999)), 2),
                        rectangles = matrix(rep(c(0.05, (y.scaling[2]-0.25)/1000*10), 10), ncol = 2, nrow = 10, byrow = TRUE),
                        lty = 0, 
                        inches = FALSE,
                        bg = "white",
                        add = TRUE)
      #
      graphics::text(x = x_positions + 0.15,
                     y = sapply(gradient_text_labels, function(x) return(0.25+round((as.numeric(x)-min(numerical_values))/abs(min(numerical_values) - max(numerical_values))*1000)*(abs(y.scaling[2])-0.5)/999)),
                     labels = gradient_text_labels,
                     adj = 0,
                     cex = 0.7)
    }
  }
  
  
  # 9. Add main and sub titles to the plot
  
  #
  if(main.title != ""){graphics::text(x = 0,
                                      y = max(y.scaling)+0.25,
                                      adj = 0.5,
                                      labels = main.title,
                                      font = 2)}
  
  #
  if(sub.title != ""){graphics::text(x = 0,
                                     y = max(y.scaling)+0.125,
                                     adj = 0.5,
                                     labels = sub.title,
                                     font = 3,
                                     cex = 0.75)}
  
  
  # 10. 
}
