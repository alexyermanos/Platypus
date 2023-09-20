#' Quantify the distance to the germline node for different features based on a Platypus Antibodyforests object.


#' @description Function to analyse and quantify the data available in an AntibodyForests object. The user can specify whether a certain number of clonotypes should be analysed. To quantify the distance either the number of nodes to germline or the
#' levenshtein distance can be used. Graphs with boxplots for node distance/Levenshtein distance to germline (including all of the nodes) and mean/median node distance/Levenshtein distance per clonotype to germline will be generated. The plots for mean/median
#' node distance and mean/median Levenshtein distance will only display the features given in the feature.ids parameter, while for the node distance/Levenshtein distance per node to germline all of the unique features that are present in the feature.id.column are
#' shown. In addition, a graph can be made in which features are specified as two groups, which shows the difference in either mean node distance to germline per clonotype or mean levenshtein distance per clonotype to germline between these two groups.The function
#' includes some parameters to customize the output, f.e. parameters to only quantify trees that contain specific features.
#'
#' @param AntibodyForest.obj Antibodyforests object as obtained from the Antibodyforests function in Platypus.
#' @param sample.ids vector containing the components of the sample_id column in your VDJ object.
#' @param feature.id.column string denoting the column in the node features data frame of the AntibodyForests object (each clonotype has it's own node features data frame, but the column name stays the same) in which the feature ids you want to analyse are given.
#' @param feature.ids vector containing the feature ids, f.e. c("s1", "s2", "s3", "s4", "s5", "s6") or f.e. c("Spleen cells", "IgM secreting BM plasma cells", "IgA secreting BM plasma cells", "GCs", "Differentiating GCs"), depending on what is present in the feature.id.column you want to analyse.
#' For the features that are given here, graphs will be generated for mean node distance and mean Levenshtein distance per clonotype per feature id.
#' @param select.specific.trees boolean denoting whether you want to select specific lineage trees, this can be used if you want to f.e. look at all trees containing "s1". If TRUE, this parameter has to be used in combination with the parameter trees.with.feature.ids.
#' @param trees.with.feature.ids vector specifying which lineage trees you want to look at, if you want to quantify all trees that f.e. contain "s1" or "s2" , you can fill in c("s1", "s2").
#' @param contain.all.feature.ids boolean specifying whether you want the trees that are used to contain all of the feature ids specified in trees.with.feature.ids. Default = FALSE, in this case all trees containing any of the features specified in trees.with.feature.ids. will be used. if TRUE,
#' only trees that contain all of the features specified in trees.with.feature.ids will be used. If you want to select trees that only contain the features specified in trees.with.feature.ids, use the parameter only.contain.feature.ids.
#' @param only.contain.feature.ids boolean denoting whether you only want to select trees containing all of the features specified in the trees.with.feature.ids parameter. Default = FALSE.
#' @param mean.or.median string denoting whether you want the function to calculate the means or the medians of the data. Default = "mean".
#' @param specify.clones boolean, is used to say whether the user wants to specify a number of clones or not. Default = FALSE, in which cases all clonotypes will be included in the function. If TRUE, the num.clones parameter can be used to specify a number of clones.
#' @param num.clones vector containing the clones which have to be analysed, f.e.if  c(1, 200) is filled in clonotype 1 until clonotype 200 will be analysed.
#' @param combine.features.line.graph boolean denoting whether the user wants to add graphs in which the features are combined into groups to show the differences between groups. Default = FALSE. If TRUE, the group1.to.combine and group2.to.combine parameters have to be used to specify the groups.
#' @param group1.to.combine vector denoting the sample ids that have to be combined to form group 1. F.e. c("s1", "s2") will combine s1 and s2 into a group. Can be used if combine.features.line.graph = TRUE.
#' @param group2.to.combine vector denoting the sample ids that have to be combined to form group 2. F.e. c("s3", "s4") will combine s3 and s4 into a group. Can be used if combine.features.line.graph = TRUE.
#'
#'@return list with 5 or 7 graphs:
#'        1) Barplot graph showing the number of nodes per feature for every feature present in the feature.id.column (only including the selected trees if select.specific.trees = TRUE)
#'        2) boxplot graph showing the node distance to germline per node (each dot is a node) for every feature present in the feature.id.column
#'        3) boxplot graph showing the mean/median node distance to germline per clonotype (each dot is a clonotype) per feature. Only features selected in the feature.ids parameter will be shown.
#'        4) boxplot graph showing the Levenshtein distance to gefrmline per node (each dot is a node) for every feature present in the feature.id.column
#'        5) boxplot graph showing the mean/median Levenshtein distance to germline per clonotype (each dot is a clonotype) per feature. Only features selected in the feature.ids parameter will be shown.
#'        In case combine.features.line.graph = TRUE:
#'        6) boxplot graph showing the mean/median node distance to germline per clonotype (each dot is a clonotype) for the two groups that were selected in the parameters group1.to.combine and group2.to.combine.
#'        7) boxplot graph showing the mean/median Levenshtein distance to germline per clonotype (each dot is a clonotype) for the two groups that were selected in the parameters group1.to.combine and group2.to.combine.
#'
#'@examples
#' \dontrun{
#' example1 <- quantify_dist_to_germ(AntibodyForest.obj = AF_VDJ_VJ_full_clust,
#' sample.ids = "s1", feature.ids = c("s1", "s2", "s3", "s4", "s5", "s6"),
#' feature.id.column = "organ_id",  select.specific.trees = FALSE,
#' mean.or.median = "mean", specify.clones = TRUE, num.clones = c(1, 20),
#' combine.features.line.graph = TRUE, group1.to.combine = c("s1", "s2"), group2.to.combine = c("s4"))
#' example1
#' }
#' \dontrun{
#' example2 <- AntibodyForests_quantify_dist_to_germ(AntibodyForest.obj = AF_VDJ_VJ_full_clust,
#' sample.ids = "s1", feature.id.column = "cluster_annotation",
#' feature.ids = c("Spleen cells", "IgM secreting BM plasma cells", "IgA secreting BM plasma cells",
#'                "GCs", "Differentiating GCs"), select.specific.trees = TRUE,
#' trees.with.feature.ids = c("GCs", "Differentiating GCs"), contain.all.feature.ids = TRUE,
#'                only.contain.feature.ids =  FALSE, mean.or.median = "mean",
#'                specify.clones = TRUE, num.clones = c(1, 20),
#'                combine.features.line.graph = FALSE)
#' example2
#' }


AntibodyForests_quantify_dist_to_germ <- function(AntibodyForest.obj, sample.ids, feature.id.column, feature.ids, select.specific.trees, trees.with.feature.ids, contain.all.feature.ids, only.contain.feature.ids, mean.or.median, specify.clones, num.clones, combine.features.line.graph, group1.to.combine, group2.to.combine) {

  if(missing(AntibodyForest.obj)) stop('Please input your AntibodyForest.obj')
  if(missing(sample.ids)) stop('Please input your sample.ids you want to analyse')
  if(missing(feature.id.column)) stop('Please input your feature.id.column, which contains the features you want to analyse within your sample.ids')
  if(missing(feature.ids)) stop('Please input your the feature.ids you want to include in the mean Levenshtein/node distance per clonotype from germline plots')
  if(missing(select.specific.trees)) select.specific.trees <- FALSE
  if(missing(trees.with.feature.ids) & select.specific.trees == TRUE) stop("Please specify the trees.with.feature.ids so that we can select the trees that contain the specific feature.ids you want to see")
  if(missing(contain.all.feature.ids)) contain.all.feature.ids <- FALSE
  if(missing(only.contain.feature.ids)) only.contain.feature.ids <- FALSE
  if(missing(mean.or.median)) mean.or.median <- "mean"
  if(missing(specify.clones)) specify.clones <- FALSE
  if(missing(num.clones) & specify.clones == TRUE) stop('Please specify the number of clones you want to analyse with the num.clones parameter')
  if(missing(combine.features.line.graph)) combine.features.line.graph <- FALSE
  if(combine.features.line.graph == FALSE) group1.to.combine <- FALSE
  if(combine.features.line.graph == FALSE) group2.to.combine <- FALSE

  feature_id <- NULL
  nodes_to_germline <- NULL
  lv_to_germline <- NULL
  value <- NULL
  group <- NULL
  clonotype_id <- NULL

  #make nested list to add all graphs per sample in
  sample.graphs.combined <- list()

  for (sample in sample.ids) {

    sample.id <- sample

    ####### Number of nodes to germline ##################################################################################

    #matrix to store number of cells per feature in
    number.of.nodes.pf <- matrix(data = 0, nrow = 1, ncol = length(feature.ids))
    colnames(number.of.nodes.pf) <- feature.ids

    #matrix to store the node distance to germline (per node) and feature info in
    nodes.to.germ.matrix <- matrix(data = NA, nrow = 0, ncol = 3)

    # make a matrix in which we can store the mean node distance to germline per feature per clonotype
    mean.nodes.to.germ.pc.matrix <- matrix(data = NA, nrow = length(AntibodyForest.obj[[sample.id]]), ncol = length(feature.ids))
    rownames(mean.nodes.to.germ.pc.matrix) <-  names(AntibodyForest.obj[[sample.id]])
    colnames(mean.nodes.to.germ.pc.matrix) <- feature.ids

    # make a matrix in which we can store the median node distance to germline per feature per clonotype
    median.nodes.to.germ.pc.matrix <- matrix(data = NA, nrow = length(AntibodyForest.obj[[sample.id]]), ncol = length(feature.ids))
    rownames(median.nodes.to.germ.pc.matrix) <-  names(AntibodyForest.obj[[sample.id]])
    colnames(median.nodes.to.germ.pc.matrix) <- feature.ids


    ####### Levenshtein distance ##########################################################################################################

    #matrix to store the levenshtein distance to germline (per node) and feature info in
    lv.to.germ.matrix <- matrix(data = NA, nrow = 0, ncol = 3)

    #matrix to store mean per feature.id per clonotype
    mean.lv.to.germ.pc.matrix <- matrix(data = NA, nrow = length(AntibodyForest.obj[[sample.id]]), ncol = length(feature.ids))
    rownames(mean.lv.to.germ.pc.matrix) <-  names(AntibodyForest.obj[[sample.id]])
    colnames(mean.lv.to.germ.pc.matrix) <- feature.ids

    #matrix to store median per feature.id per clonotype
    median.lv.to.germ.pc.matrix <- matrix(data = NA, nrow = length(AntibodyForest.obj[[sample.id]]), ncol = length(feature.ids))
    rownames(median.lv.to.germ.pc.matrix) <-  names(AntibodyForest.obj[[sample.id]])
    colnames(median.lv.to.germ.pc.matrix) <- feature.ids


    ####### line graph ##########################################################################################################

    # matrix for line graph part - depending on combine.features.line graph it will either have two or three columns, for nodes and for levenshtein distance
    lg.nodes.groups.combined.matrix <- matrix(data = NA, nrow = 0, ncol = 3)
    colnames(lg.nodes.groups.combined.matrix) <- c("clonotype_id", "group1", "group2")

    lg.lv.groups.combined.matrix <- matrix(data = NA, nrow = 0, ncol = 3)


    ####### specify.clones - num.clones ##########################################################################################################

    # making sure that the right number of clones goes through the function
    total.clones <- as.numeric(length(AntibodyForest.obj[[sample.id]]))


    if(specify.clones == TRUE){
      start.clone <- num.clones[1]
      stop.clone <- num.clones[2]

      if (stop.clone>total.clones){
        message("Wanted clones exceed the clones in the VDJ. Computing ref for clone ", num.clones[1], " until clone ", total.clones, " instead")
        start.clone = num.clones[1]
        stop.clone = total.clones
      }

    }else if (specify.clones == FALSE){
      message("Computing ref for all ", total.clones," clones")
      start.clone <- 1
      stop.clone <- total.clones
    }

    #################################################################################################################################################
    #For loop #######################################################################################################################################


    # make a for loop that goes over the clonotypes in the AntibodyForest Object

    for (clonotype.num in start.clone : stop.clone) {

      clonotype.id <- paste("clonotype", clonotype.num, sep = "")

      # Access the clonotype using the identifier
      if (clonotype.id %in% names(AntibodyForest.obj[[sample.id]])) {
        clonotype.info <- AntibodyForest.obj[[sample.id]][[clonotype.id]]

        # Extract node features to clonotype
        node.features <- clonotype.info@node_features

        # If the clonotype is not present in the antibody forest object, the function will skip over that one.
      } else {
        message("Clonotype", clonotype.id, " not found. Skipping to the next clonotype.")
        next
      }

      #if select.specific.trees == TRUE, there are three different options.
      if (select.specific.trees == TRUE) {

        unique.node.features <- as.character(unique(node.features[, feature.id.column]))
        unique.node.features <- unique.node.features[unique.node.features != "germline"]


        if (contain.all.feature.ids == TRUE) {

          if (only.contain.feature.ids == TRUE) {


            if ((identical(sort(unique.node.features), sort(trees.with.feature.ids)))) {

              message("Clonotype", clonotype.id, " consists of all the features in the trees.with.feature.ids parameter and will be used")

            } else {

              message("Clonotype", clonotype.id, " does not consist only of the features in the feature.ids parameter. Skipping to the next clonotype.")
              next
            }

          } else if (only.contain.feature.ids == FALSE) {

            if (all(trees.with.feature.ids %in% sort(unique.node.features))) {

              message("Clonotype", clonotype.id, " contains all the features in the trees.with.feature.ids parameter and will be used")

            } else {

              message("Clonotype", clonotype.id, " does not contain all of the features in the trees.with.feature.ids parameter. Skipping to the next clonotype.")

              next

            }

          }

        } else if (contain.all.feature.ids == FALSE) {


          if (sum(unique.node.features %in% trees.with.feature.ids) > 1) {

            message("Clonotype", clonotype.id, " contains one or multiple features that are present in the trees.with.feature.ids parameter and will be used")

          } else if (sum(unique.node.features %in% trees.with.feature.ids) == 0) {

            message("Clonotype", clonotype.id, " does not contain any of the features in the trees.with.feature.ids parameter. Skipping to the next clonotype.")
            next

          }


        }


      }


      ####### Number of nodes per feature ##################################################################################


      for (i in feature.ids) {

        right.column <- node.features[, feature.id.column]
        number.of.nodes.pf[1 , i] <- number.of.nodes.pf[1 , i] + length(right.column[right.column == i])

      }



      ####### Number of nodes to germline ##################################################################################

      #add a column to node_features to add the number of nodes to germline
      node.features$nodes_to_germline <- NA

      for (j in node.features$label) {

        #calculate distance to the germline (which is the last node)
        node.dist.to.germ <- igraph::distances(clonotype.info@tree, v = j, to = utils::tail(node.features$label, 1), weights = NA, algorithm = "unweighted")

        node.features[node.features$label == j, "nodes_from_germline"] <- node.dist.to.germ

      }


      #extract feature id column and nodes to germline column to node.features
      feature.id.pc <- node.features[, feature.id.column]
      nodes.to.germline.pc <- node.features$nodes_from_germline

      #bind the feature id column and nodes to germline column with the clonotype.id and add colnames
      nodes.to.germ.pc.matrix <- cbind(feature.id.pc, nodes.to.germline.pc, clonotype.id)
      colnames(nodes.to.germ.pc.matrix) <- c("feature_id", "nodes_to_germline", "clonotype_id")

      #add rows to nodes.to.germ.matrix
      nodes.to.germ.matrix <- rbind(nodes.to.germ.matrix, nodes.to.germ.pc.matrix)

      ####### mean per feature ##########
      if (mean.or.median == "mean") {

        #make data.frame because matrix does not work for some reason for mean part
        nodes.to.germ.pc.df <- as.data.frame(nodes.to.germ.pc.matrix)

        #make sure that function recognizes numbers in nodes.to.germ.pc.df$nodes_to_germline as numbers
        nodes.to.germ.pc.df$nodes_to_germline <- as.numeric(nodes.to.germ.pc.df$nodes_to_germline)

        #make a matrix to store the mean levenshtein distance per feature per clonotype in
        mean.nodes.to.germ.pf.pc <- matrix(NA, nrow = 1, ncol = length(feature.ids))
        colnames(mean.nodes.to.germ.pf.pc) <- feature.ids

        # loop over the features
        for (j in feature.ids) {

          mean <- mean(nodes.to.germ.pc.df[nodes.to.germ.pc.df$feature_id == j, "nodes_to_germline"])

          mean.nodes.to.germ.pf.pc[1, j] <- mean

        }

        #add the mean of the faetures per clonotype to the matrix with all the means of the features
        mean.nodes.to.germ.pc.matrix[clonotype.id, ] <- mean.nodes.to.germ.pf.pc

        ####### median per feature ##########
      } else if (mean.or.median == "median") {

        #make data.frame because matrix does not work for some reason for median part
        nodes.to.germ.pc.df <- as.data.frame(nodes.to.germ.pc.matrix)

        #make sure that function recognizes numbers in nodes.to.germ.pc.df$nodes_to_germline as numbers
        nodes.to.germ.pc.df$nodes_to_germline <- as.numeric(nodes.to.germ.pc.df$nodes_to_germline)

        #make a matrix to store the mean levenshtein distance per feature per clonotype in
        median.nodes.to.germ.pf.pc <- matrix(NA, nrow = 1, ncol = length(feature.ids))
        colnames(median.nodes.to.germ.pf.pc) <- feature.ids

        # loop over the features
        for (j in feature.ids) {

          median <- median(nodes.to.germ.pc.df[nodes.to.germ.pc.df$feature_id == j, "nodes_to_germline"])

          median.nodes.to.germ.pf.pc[1, j] <- median

        }

        #add the mean of the faetures per clonotype to the matrix with all the means of the features
        median.nodes.to.germ.pc.matrix[clonotype.id, ] <- median.nodes.to.germ.pf.pc


      }


      ####### Levenshtein distance ##########################################################################################################

      #take distance_to_germline and feature_id column out of the node.features data frame and cbind together with the clonotype.id
      lv.dist.to.germ <- node.features[, "distance_from_germline"]
      feature.info <- node.features[, feature.id.column]
      lv.to.germ.pc <- cbind(lv.dist.to.germ, feature.info, clonotype.id)

      #add data per clonotype to lv.to.germ.matrix
      lv.to.germ.matrix <- rbind(lv.to.germ.matrix, lv.to.germ.pc)


      ####### mean per feature ##########
      if (mean.or.median == "mean") {

        #make a matrix to store the mean Levenshtein distance per feature per clonotype in
        mean.lv.to.germ.pf.pc <- matrix(NA, nrow = 1, ncol = length(feature.ids))
        colnames(mean.lv.to.germ.pf.pc) <- feature.ids

        # loop over the features to get the mean Levenshtein distance per feature per clonotype
        for (j in feature.ids) {

          mean <- mean(node.features[node.features[,feature.id.column] == j, "distance_from_germline"])

          mean.lv.to.germ.pf.pc[1, j] <- mean

        }

        #add the mean of the features per clonotype to the matrix with all the means of the features
        mean.lv.to.germ.pc.matrix[clonotype.id, ] <- mean.lv.to.germ.pf.pc

      } else if (mean.or.median == "median") {

        #make a matrix to store the mean Levenshtein distance per feature per clonotype in
        median.lv.to.germ.pf.pc <- matrix(NA, nrow = 1, ncol = length(feature.ids))
        colnames(median.lv.to.germ.pf.pc) <- feature.ids

        # loop over the features to get the mean Levenshtein distance per feature per clonotype
        for (j in feature.ids) {

          median <- median(node.features[node.features[,feature.id.column] == j, "distance_from_germline"])

          median.lv.to.germ.pf.pc[1, j] <- median

        }

        #add the mean of the features per clonotype to the matrix with all the means of the features
        median.lv.to.germ.pc.matrix[clonotype.id, ] <- median.lv.to.germ.pf.pc

      }


      ####### line graph ##########################################################################################################

      ## first select if the groups have to be combined or not
      if (combine.features.line.graph == TRUE) {

        ##### Node distance ##########################

        #first make a data frame to the nodes.to.germ.pc.matrix that I made in the Number of nodes to germline part
        nodes.to.germ.pc.df <- as.data.frame(nodes.to.germ.pc.matrix)
        colnames(nodes.to.germ.pc.df) <- c("feature_id", "nodes_to_germline", "clonotype_id")
        nodes.to.germ.pc.df$nodes_to_germline <- as.numeric(nodes.to.germ.pc.df$nodes_to_germline)

        #filter for group 1 and for group 2
        nodes.group.1.df <- nodes.to.germ.pc.df %>%
          dplyr::filter(feature_id %in% group1.to.combine)

        nodes.group.2.df <- nodes.to.germ.pc.df %>%
          dplyr::filter(feature_id %in% group2.to.combine)

        if (mean.or.median == "mean") {

          #calculate mean for group 1
          nodes.group1 <- nodes.group.1.df %>%
            dplyr::summarise(mean_value = mean(nodes_to_germline))

          #calculate mean group 2
          nodes.group2 <- nodes.group.2.df %>%
            dplyr::summarise(mean_value = mean(nodes_to_germline))

        } else if (mean.or.median == "median") {

          #calculate median for group 1
          nodes.group1 <- nodes.group.1.df %>%
            dplyr::summarise(median_value = median(nodes_to_germline))

          #calculate mean group 2
          nodes.group2 <- nodes.group.2.df %>%
            dplyr::summarise(median_value = median(nodes_to_germline))

        }

        #cbind values with clonotype.id
        row.nodes.groups.combined.pc.matrix <- cbind(clonotype.id, as.numeric(nodes.group1), as.numeric(nodes.group2))

        #add data per clonotype to lg.nodes.groups.combined.matrix
        lg.nodes.groups.combined.matrix <- rbind(lg.nodes.groups.combined.matrix, row.nodes.groups.combined.pc.matrix )


        ##### Levenshtein distance #####################

        lv.to.germ.pc.df <- as.data.frame(lv.to.germ.pc)
        colnames(lv.to.germ.pc.df) <- c("lv_to_germline", "feature_id", "clonotype_id")
        lv.to.germ.pc.df$lv_to_germline <- as.numeric(lv.to.germ.pc.df$lv_to_germline)

        #filter for group 1 and for group 2
        lv.group.1.df <- lv.to.germ.pc.df %>%
          dplyr::filter(feature_id %in% group1.to.combine)

        lv.group.2.df <- lv.to.germ.pc.df %>%
          dplyr::filter(feature_id %in% group2.to.combine)

        if (mean.or.median == "mean") {

          #calculate mean for group 1
          lv.group1 <- lv.group.1.df %>%
            dplyr::summarise(mean_value = mean(lv_to_germline))

          #calculate mean group 2
          lv.group2 <- lv.group.2.df %>%
            dplyr::summarise(mean_value = mean(lv_to_germline))

        } else if (mean.or.median == "median") {

          #calculate median for group 1
          lv.group1 <- lv.group.1.df %>%
            dplyr::summarise(median_value = median(lv_to_germline))

          #calculate mean group 2
          lv.group2 <- lv.group.2.df %>%
            dplyr::summarise(median_value = median(lv_to_germline))
        }

        row.lv.groups.combined.pc.matrix <- cbind(clonotype.id, lv.group1, lv.group2)
        colnames(row.lv.groups.combined.pc.matrix) <- c("clonotype_id", "group1", "group2")

        #add data to lg.lv.groups.combined.matrix
        lg.lv.groups.combined.matrix <- rbind(lg.lv.groups.combined.matrix, row.lv.groups.combined.pc.matrix )




      }

    }

    #####################################################################################################################################################
    #End For loop #######################################################################################################################################


    ###### Number of cells per feature plot ################################################################################################################

    nr.of.nodes.pf.df <- as.data.frame(number.of.nodes.pf)

    nr.of.nodes.pf.df.long <- tidyr::gather(nr.of.nodes.pf.df, key = "feature_id", value = "value", 1:length(feature.ids))

    if (specify.clones == TRUE) {
      title.nr.of.nodes.pf.plot <- paste("Number of nodes per feature for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
    } else if (specify.clones == FALSE) {
      title.nr.of.nodes.pf.plot <- paste("Number of nodes per feature for", sample.id)
    }

    nr.of.nodes.pf.plot <- ggplot2::ggplot(data = nr.of.nodes.pf.df.long, ggplot2::aes(x= feature_id, y = value)) +
      ggplot2::geom_bar(stat="identity", fill= "seagreen") +
      ggplot2::geom_text(ggplot2::aes(label=value), vjust=-0.3, color="black", size=7)+
      ggplot2::theme_minimal() +
      ggplot2::labs(title = title.nr.of.nodes.pf.plot, x = " ", y = "Number of nodes")


    ####### Node distance plots ################################################################################################################

    #first let's make a  data frame with only the feature_ids and the nodes_to_germline column from the earlier made matrix
    nodes.to.germ.df <- as.data.frame(nodes.to.germ.matrix)
    colnames(nodes.to.germ.df) <- c("feature_id", "nodes_to_germline", "clonotype_id")

    #R does not recognize the numbers as numbers and the characters as characters
    nodes.to.germ.df$feature_id <- as.character(nodes.to.germ.df$feature_id)
    nodes.to.germ.df$nodes_to_germline <- as.numeric(nodes.to.germ.df$nodes_to_germline)

    if (specify.clones == TRUE) {
      title.bp.nodes.to.germ <- paste("Node distance to germline per node for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
    } else if (specify.clones == FALSE) {
      title.bp.nodes.to.germ <- paste("Node distance to germline per node for", sample.id)
    }

    boxplot.nodes.to.germline <- ggplot2::ggplot(data = nodes.to.germ.df, ggplot2::aes(x=feature_id, y=nodes_to_germline, fill=feature_id)) +
      ggplot2::geom_boxplot()  +
      ggplot2::geom_jitter(color="black", size=0.4, alpha=0.9) +
      ggplot2::scale_fill_brewer(palette="PiYG") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10)) +
      ggplot2::labs(title = title.bp.nodes.to.germ, x = " ", y = "Node distance to germline")


    ####### mean or median per feature ##########

    #specify whether we want to visualize the mean or median
    if (mean.or.median == "mean") {

      nodes.to.germ.pc.df <- as.data.frame(mean.nodes.to.germ.pc.matrix)

    } else if (mean.or.median == "median") {

      nodes.to.germ.pc.df <- as.data.frame(median.nodes.to.germ.pc.matrix)

    }

    #create extra column with clonotype_ids
    nodes.to.germ.pc.df$clonotype_id <- rownames(nodes.to.germ.pc.df)

    colnames.to.gather.nodes <- colnames(nodes.to.germ.pc.df)[-length(colnames(nodes.to.germ.pc.df))]
    nodes.to.germ.pc.df.long <- tidyr::gather(nodes.to.germ.pc.df, key = "feature_id", value = "value", colnames.to.gather.nodes)
    nodes.to.germ.pc.df.long <- stats::na.omit(nodes.to.germ.pc.df.long)

    if (mean.or.median == "mean") {

      if (specify.clones == TRUE) {
        title.bp.nodes.to.germ <- paste("Mean node distance to germline per clonotype for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
      } else if (specify.clones == FALSE) {
        title.bp.nodes.to.germ <- paste("Mean node distance to germline per clonotype for", sample.id)
      }

    } else if (mean.or.median == "median") {

      if (specify.clones == TRUE) {
        title.bp.nodes.to.germ <- paste("Median node distance to germline per clonotype for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
      } else if (specify.clones == FALSE) {
        title.bp.nodes.to.germ <- paste("Median node distance to germline per clonotype for", sample.id)
      }

    }

    boxplot.mean.or.median.nodes.pc <- ggplot2::ggplot(data = nodes.to.germ.pc.df.long, ggplot2::aes(x=feature_id, y=value, fill=feature_id)) +
      ggplot2::geom_boxplot()  +
      ggplot2::geom_jitter(color="black", size=0.4, alpha=0.9) +
      ggplot2::scale_fill_brewer(palette="Paired") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = title.bp.nodes.to.germ, x = " ", y = "Node distance to germline")


    ######### Levenshtein distance plots ################################################################################################################################

    lv.to.germ.df <- as.data.frame(lv.to.germ.matrix)
    colnames(lv.to.germ.df) <- c("lv_to_germline", "feature_id", "clonotype_id")

    lv.to.germ.df$feature_id <- as.character(lv.to.germ.df$feature_id)
    lv.to.germ.df$lv_to_germline <- as.numeric(lv.to.germ.df$lv_to_germline)

    if (specify.clones == TRUE) {
      title.bp.lv.to.germ <- paste("Levenshtein distance to germline per node for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
    } else if (specify.clones == FALSE) {
      title.bp.lv.to.germ <- paste("Levenshtein distance to germline per node for", sample.id)
    }

    boxplot.lv.to.germline <- ggplot2::ggplot(data = lv.to.germ.df, ggplot2::aes(x=feature_id, y=lv_to_germline, fill=feature_id)) +
      ggplot2::geom_boxplot()  +
      ggplot2::geom_jitter(color="black", size=0.4, alpha=0.9) +
      ggplot2::scale_fill_brewer(palette="PiYG") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10)) +
      ggplot2::labs(title = title.bp.lv.to.germ, x = " ", y = "Levenshtein distance to germline")

    ####### mean or median per feature ##########

    if (mean.or.median == "mean") {

      lv.to.germ.pc.df <- as.data.frame(mean.lv.to.germ.pc.matrix)

    } else if (mean.or.median == "median") {

      lv.to.germ.pc.df <- as.data.frame(median.lv.to.germ.pc.matrix)

    }

    #add clonotype_ids as a column with rows instead of just rownames
    lv.to.germ.pc.df$clonotype_id <- rownames(lv.to.germ.pc.df)

    colnames.to.gather.lv <- colnames(lv.to.germ.pc.df)[-length(colnames(lv.to.germ.pc.df))]
    lv.to.germ.pc.df.long <- tidyr::gather(lv.to.germ.pc.df, key = "feature_id", value = "value", colnames.to.gather.lv)
    lv.to.germ.pc.df.long <- stats::na.omit(lv.to.germ.pc.df.long)


    #make sure right title is added to the boxplot
    if (mean.or.median == "mean") {

      if (specify.clones == TRUE) {
        title.bp.lv.to.germ <- paste("Mean Levenshtein distance to germline per clonotype for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
      } else if (specify.clones == FALSE) {
        title.bp.lv.to.germ <- paste("Mean Levenshtein distance to germline per clonotype for", sample.id)
      }

    } else if (mean.or.median == "median") {

      if (specify.clones == TRUE) {
        title.bp.lv.to.germ <- paste("Median Levenshtein distance to germline per clonotype for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
      } else if (specify.clones == FALSE) {
        title.bp.lv.to.germ <- paste("Median Levenshtein distance to germline per clonotype for", sample.id)
      }

    }

    boxplot.mean.or.median.lv.pc <- ggplot2::ggplot(data = lv.to.germ.pc.df.long, ggplot2::aes(x=feature_id, y=value, fill=feature_id)) +
      ggplot2::geom_boxplot()  +
      ggplot2::geom_jitter(color="black", size=0.4, alpha=0.9) +
      ggplot2::scale_fill_brewer(palette="Paired") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = title.bp.lv.to.germ, x = " ", y = "Levenshtein distance to germline")



    ######### Line graph ################################################################################################################################

    ## first select if the groups have to be combined or not
    if (combine.features.line.graph == TRUE) {

      ##### Node distance ########

      #make a long data.frame to lg.groups.combined.matrix
      nodes.groups.combined.df <- as.data.frame(lg.nodes.groups.combined.matrix)
      nodes.groups.combined.df$group1 <- as.numeric(nodes.groups.combined.df$group1)
      nodes.groups.combined.df$group2 <- as.numeric(nodes.groups.combined.df$group2)

      if(all(is.na(nodes.groups.combined.df$group1)) & all(is.na(nodes.groups.combined.df$group2))) {

        lg.nodes.boxplot.groups.combined <- "There are no values available for these groups in this sample"

      } else {
        #take out rows with NA
        nodes.groups.combined.df <- stats::na.omit(nodes.groups.combined.df)


        nodes.groups.combined.df.long <- nodes.groups.combined.df %>%
          tidyr::pivot_longer(cols = tidyselect::starts_with("group"),  # Select columns starting with "means_group"
                       names_to = "group",                # New column for the group names
                       values_to = "value")

        #make sure it recognizes the means as numerical values
        nodes.groups.combined.df.long$value <- as.numeric(nodes.groups.combined.df.long$value)

        #make labels to be showed on x axis
        group1.name <- paste(group1.to.combine, collapse = ", ")
        group2.name <- paste(group2.to.combine, collapse = ", ")


        #show the right title on the graph
        if (mean.or.median == "mean") {

          if (specify.clones == TRUE) {
            title.lg <- paste("Mean node distance to germline per group for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
          } else if (specify.clones == FALSE) {
            title.lg <- paste("Mean node distance to germline per group for", sample.id)
          }
          ylab <- "Mean node distance to germline"

        } else if (mean.or.median == "median") {

          if (specify.clones == TRUE) {
            title.lg <- paste("Median node distance to germline per group for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
          } else if (specify.clones == FALSE) {
            title.lg <- paste("Median node distance to germline per group for", sample.id)
          }
          ylab <- "Median node distance to germline"

        }


        lg.nodes.boxplot.groups.combined <- ggplot2::ggplot(data = nodes.groups.combined.df.long, ggplot2::aes(x=group, y=value, fill=group)) +
          ggplot2::geom_boxplot()  +
          ggplot2::geom_jitter(color="black", size=0.4, alpha=0.9) +
          ggplot2::geom_line(ggplot2::aes(group = clonotype_id), color = "gray", alpha = 0.5) +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = title.lg , x = " ", y = ylab) +
          ggplot2::scale_x_discrete(labels = c(group1.name, group2.name))
      }


      ##### Levenshtein distance ########

      lg.lv.groups.combined.df <- as.data.frame(lg.lv.groups.combined.matrix)

      if(all(is.na(lg.lv.groups.combined.df$group2)) & all(is.na(lg.lv.groups.combined.df$group2))) {

        lg.lv.boxplot.groups.combined <- "There are no values available for these groups in this sample"

      } else {

        #take out rows with NA
        lg.lv.groups.combined.df <- stats::na.omit(lg.lv.groups.combined.df)

        #make the data frame long
        lv.groups.combined.df.long <- lg.lv.groups.combined.df %>%
          tidyr::pivot_longer(cols = tidyselect::starts_with("group"),  # Select columns starting with "means_group"
                       names_to = "group",                # New column for the group names
                       values_to = "value")

        #make sure it recognizes the means as numerical values
        lv.groups.combined.df.long$value <- as.numeric(lv.groups.combined.df.long$value)

        if (mean.or.median == "mean") {

          if (specify.clones == TRUE) {
            title.lg <- paste("Mean Levenshtein distance to germline per group for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
          } else if (specify.clones == FALSE) {
            title.lg <- paste("Mean Levenshtein distance to germline per group for", sample.id)
          }

          ylab <- "Mean Levenshtein distance to germline"

        } else if (mean.or.median == "median") {

          if (specify.clones == TRUE) {
            title.lg <- paste("Median Levenshtein distance to germline per group for clonotype", start.clone, "till clonotype", stop.clone, "for", sample.id)
          } else if (specify.clones == FALSE) {
            title.lg <- paste("Median Levenshtein distance to germline per group for", sample.id)
          }

          ylab <- "Median Levenshtein distance to germline"

        }

        group1.name <- paste(group1.to.combine, collapse = ", ")
        group2.name <- paste(group2.to.combine, collapse = ", ")

        lg.lv.boxplot.groups.combined <- ggplot2::ggplot(data = lv.groups.combined.df.long, ggplot2::aes(x=group, y=value, fill=group)) +
          ggplot2::geom_boxplot()  +
          ggplot2::geom_jitter(color="black", size=0.4, alpha=0.9) +
          ggplot2::geom_line(ggplot2::aes(group = clonotype_id), color = "gray", alpha = 0.5) +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = title.lg, x = " ", y = ylab) +
          ggplot2::scale_x_discrete(labels = c(group1.name, group2.name))
      }

    }



    ############# Return the graphs #####################################################################################################################

    list.graphs.pf <- list()
    if (combine.features.line.graph == TRUE) {

      list.graphs.pf$nr.of.nodes.pf.plot <- nr.of.nodes.pf.plot
      list.graphs.pf$boxplot.nodes.to.germline <- boxplot.nodes.to.germline
      list.graphs.pf$boxplot.mean.or.median.nodes.pc <- boxplot.mean.or.median.nodes.pc
      list.graphs.pf$boxplot.lv.to.germline <- boxplot.lv.to.germline
      list.graphs.pf$boxplot.mean.or.median.lv.pc <- boxplot.mean.or.median.lv.pc
      list.graphs.pf$lg.nodes.boxplot.groups.combined <- lg.nodes.boxplot.groups.combined
      list.graphs.pf$lg.lv.boxplot.groups.combined <- lg.lv.boxplot.groups.combined

    } else {

      list.graphs.pf$nr.of.nodes.pf.plot <- nr.of.nodes.pf.plot
      list.graphs.pf$boxplot.nodes.to.germline <- boxplot.nodes.to.germline
      list.graphs.pf$boxplot.mean.nodes.pc <- boxplot.mean.or.median.nodes.pc
      list.graphs.pf$boxplot.lv.to.germline <- boxplot.lv.to.germline
      list.graphs.pf$boxplot.mean.lv.pc <- boxplot.mean.or.median.lv.pc


    }


    # Add the inner list to the outer list
    sample.graphs.combined[[sample.id]] <- list.graphs.pf

  }

  return(sample.graphs.combined)
}


