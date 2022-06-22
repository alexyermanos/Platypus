
AntibodyForests_label_propagation <- function(trees,
                                              features,
                                              propagation.algorithm,
                                              diffusion.n.iter,
                                              diffusion.threshold,
                                              parallel){


  if(missing(trees)) stop('Please input a nested list of AntibodyForests objects.')
  if(missing(features)) features <- NULL
  if(missing(propagation.algorithm)) propagation.algorithm <- 'diffusion'
  if(missing(diffusion.n.iter)) diffusion.n.iter <- 20
  if(missing(diffusion.threshold)) diffusion.threshold <- 1e-3
  if(missing(parallel)) parallel <- T



  get_label_dict <- function(trees, features){
    classes_dict <- vector(mode = 'list', length = length(features))

    for(i in 1:length(features)){

      unique_values <- c()

      if(inherits(trees, 'list')){

        for(j in 1:length(trees)){

          for(k in 1:length(trees[[j]])){
            g <- trees[[j]][[k]]@tree

            vals <- unlist(igraph::vertex_attr(g, name = features[i]))

            unique_values <- c(unique_values, unique(vals))
          }
        }
      }else{
        g <- trees@tree
        vals <- unlist(igraph::vertex_attr(g, name = features[i]))
        unique_values <- unique(vals)
      }

      unique_values <- unique(unique_values)
      unique_values <- unique_values[unique_values != 'germline' & unique_values != 'unknown' & unique_values != 'intermediate']
      ids <- 1:length(unique_values)
      names(ids) <- unique_values
      classes_dict[[i]] <-  ids
    }

    names(classes_dict) <- features
    return(classes_dict)
  }


  get_feature_vector <- function(tree, label_dict, feature){

    id_dict <- label_dict[[feature]]

    g <- tree@tree

    labels <- igraph::vertex_attr(g, name = feature)

    if(paste0(feature, '_counts') %in% igraph::vertex_attr_names(g)){
      counts <- igraph::vertex_attr(g, name=paste0(feature, '_counts'))
      max_indices <- lapply(counts, function(x) which.max(x))
      labels <- unlist(mapply(function(x,y) x[y], labels, max_indices))
    }

    labels <- unlist(labels)
    labels[labels == 'intermediate' | labels == 'germline' | is.na(labels) | labels == ''] <-  'unknown'
    labels <- unname(id_dict[labels])
    labels[is.na(labels)] <- -1


    return(labels)
  }

  one_hot_features <- function(feature.vector, n.classes){

    one_hot_matrix <- matrix(0, length(feature.vector), n.classes)
    ind <- which(feature.vector != -1)


    for(i in ind){
      one_hot_matrix[i, feature.vector[i]] <- 1
    }

    return(one_hot_matrix)
  }

  label_propagation_diffusion <- function(tree, features, label_dict){

    g <- tree@tree
    g <- igraph::as.undirected(g)
    igraph::E(g)$weight <- 1/igraph::E(g)$weight
    A <- igraph::as_adjacency_matrix(g, attr = 'weight')
    D <- matrix(0, length(igraph::V(g)), length(igraph::V(g)))
    diag(D) <- igraph::degree(g)
    D_inv <- solve(D)
    T <- D_inv %*% as.matrix(A)


    for(feature in features){
      feature_vector <- tree %>% get_feature_vector(label_dict = label_dict, feature = feature)
      fixed <- feature_vector != -1
      X_init <- feature_vector %>% one_hot_features(n.classes = length(label_dict[[feature]]))
      X_1 <- X_init

      n <- 1
      current_diff <- .Machine$integer.max

      while(n <= diffusion.n.iter | current_diff > diffusion.threshold){
        X_0 <- X_1
        X_1 <- T %*% X_0
        X_1[fixed,] <- X_init[fixed,]

        current_diff <- max(abs(X_1 - X_0))
        n <- n + 1
      }

      new_labels <- unlist(lapply(1:nrow(X_1), function(x) which.max(X_1[x,])))
      dict <- label_dict[[feature]]
      new_labels <- names(dict[new_labels])
      g <- igraph::set_vertex_attr(g, name = paste0(feature, '_label_propagation'), value = new_labels)
    }

    tree@tree <- g

    return(tree)
  }


  label_propagation_communities <- function(tree, features, label_dict){

    g <- tree@tree
    igraph::E(g)$weight <- 1/igraph::E(g)$weight

    for(feature in features){
      feature_vector <- get_feature_vector(tree, label_dict = label_dict, feature = feature)
      fixed <- feature_vector != -1
      cluster <- igraph::cluster_label_prop(g, weights = NULL, initial = feature_vector, fixed = fixed)
      new_labels <- cluster$membership

      dict <- label_dict[[feature]]
      new_labels <- names(dict[new_labels])
      g <- igraph::set_vertex_attr(g, name = paste0(feature, '_label_propagation'), value = new_labels)
    }

    tree@tree <- g

    return(tree)
  }

  label_propagation <- function(tree, features, label_dict, propagation.algorithm){

    if(propagation.algorithm == 'diffusion'){
      tree <- label_propagation_diffusion(tree, features, label_dict)
    }else if(propagation.algorithm == 'communities'){
      tree <- label_propagation_communities(tree, features, label_dict)
    }else{
      stop('Label propagation algorithm not available!')
    }

    return(tree)
  }

  get_feature_names <- function(trees, features){

    if(is.null(features)){
      if(inherits(trees, 'list')){
        features <- trees[[1]][[1]]@feature_names
      }else if(inherits(trees, 'AntibodyForests')){
        features <- trees@feature_names
      }

      if(is.null(features)){
        stop('Could not find the features to perform label propagation on! Please provide the feature names in the features parameter!')
      }
    }

    return(features)
  }

  features <- get_feature_names(trees, features)
  label_dict <- get_label_dict(trees, features)

  if(inherits(trees, 'list')){

    for(i in 1:length(trees)){

      if(parallel){
        requireNamespace('parallel')
        cores <- parallel::detectCores()
        trees[[i]] <- parallel::mclapply(trees[[i]], mc.cores = cores,
                                                 FUN = function(x) {label_propagation(x, features, label_dict, propagation.algorithm)
                                                                    })

      }else{

        trees[[i]] <- lapply(trees[[i]], function(x) label_propagation(x, features, label_dict, propagation.algorithm))

      }
    }

  }else if(inherits(trees, 'AntibodyForests')){
    trees <- label_propagation(trees, features, label_dict, propagation.algorithm)


  }else{
    stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
  }

  return(trees)
}
