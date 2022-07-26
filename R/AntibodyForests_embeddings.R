#'Structural node embeddings for the AntibodyForests minimum spanning trees/ sequence similarity networks

#'@description Structural node embeddings algorithms of the AntibodyForests networks. Supported algorithms include: node2vec (https://arxiv.org/abs/1607.00653) and spectral graph embedding on either the adjacency or the Laplacian matrix.  Currently the node2vec model is supported as long as Rkeras is installed.

#' @param trees AntibodyForests object/list of AntibodyForests objects - the resulting sequence similarity or minimum spanning tree networks from the AntibodyForests function
#' @param graph.type string - the graph type available in the AntibodyForests object which will be used as the function input.
#' Currently supported network/analysis types: 'tree' (for the minimum spanning trees or sequence similarity networks obtained from the main AntibodyForests function), 'heterogeneous' for the bipartite graphs obtained via AntibodyForests_heterogeneous, 'dynamic' for the dynamic networks obtained from AntibodyForests_dynamics.
#' @param embedding.method string - the embeddings model/algorithm. 'node2vec' for an implementation of graph random walk and node2vec using R-keras (might be slow depending on graph size), 'spectral_adjacency' for spectral graph embeddings of the adjacency matrix (using igraph's embed_adjacency_matrix() function), 'spectral_laplacian' for embedding the Laplacian matrix (using igraph's embed_laplacian_matrix() function).
#' @param dim.reduction string - dimensionality reduction algorithm for the resulting node2vec embeddings. Currently implemented methods include: 'umap', 'tsne' and 'pca'.
#' @param color.by vector of strings - features to color the resulting scatter plots by. These features must be included as igraph vertex attributes when creating the AntibodyForests objects, by including them in the node.features parameter.
#' @param num.walks integer - number of biased random walks to be performed for the node2vec training dataset.
#' @param num.step integer - number of steps per biased random walk.
#' @param p numeric - probability of revisiting the same node already vistied in a random walk step (= return parameter).
#' @param q numeric - probability of 'jumping' to a node closer or farther away from the node visited at step x (e.g., q > 1, random walk is biased to closer nodes, q < 1, random walk will 'jump' to farher nodes more frequently).
#' @param window.size integer - size of sampling window in the skipgram model.
#' @param num.negative.samples integer - number of negative samples to be considered in the skipgram model.
#' @param embedding.dim integer - latent/embedding dimension of the node2vec output vectors.
#' @param batch.size integer - training batch size of the node2vec model.
#' @param epochs integer - number of training epochs for the node2vec model.
#' @param tsne.perplexity numeric - T-SNE reduction perplexity.
#' @param seed integer - random seed for the random walk steps of the node2vec model.
#' @param parallel boolean - whether to execute the random walks in parallel or not.


#' @return A scatterplot of reduced vector embeddings for each node in the graphs, colored by the features specified in color.by.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_embeddings(output_networks, graph.type = 'tree', embedding.method = 'node2vec', dim.reduction = 'pca', num.walks = 10, num.steps = 10, embedding.dim = 64, batch.size = 32, epochs = 50)
#'}


#TO DO: ADD OPTION TO SAVE THE EMBEDDING NETWORK FOR NODE2VEC


AntibodyForests_embeddings <- function(trees,
                                       graph.type,
                                       embedding.method,
                                       dim.reduction,
                                       color.by,
                                       num.walks,
                                       num.steps,
                                       p,
                                       q,
                                       window.size,
                                       num.negative.samples,
                                       embedding.dim,
                                       batch.size,
                                       epochs,
                                       tsne.perplexity,
                                       seed,
                                       parallel
                                       ){

   if(missing(trees)) stop('Please input a nested list of AntibodyForests objects or a single AntibodyForests object')
   if(missing(graph.type)) graph.type <- 'tree'
   if(missing(embedding.method)) embedding.method <- 'node2vec'
   if(missing(dim.reduction)) dim.reduction <- 'pca'
   if(missing(color.by)) color.by <- 'node_type'
   if(missing(num.walks)) num.walks <- 10
   if(missing(num.steps)) num.steps <- 10
   if(missing(p)) p <- 1
   if(missing(q)) q <- 1
   if(missing(window.size)) window.size <- 5
   if(missing(num.negative.samples)) num.negative.samples <- 4
   if(missing(embedding.dim)) embedding.dim <- 128
   if(missing(batch.size)) batch.size <- 128
   if(missing(epochs)) epochs <- 20
   if(missing(tsne.perplexity)) tsne.perplexity <- 5
   if(missing(seed)) seed <- 1
   if(missing(parallel)) parallel <- T

   rw_step <- function(g, prev_node, current_node){

     sampling_weights <- c()

     if(igraph::is_directed(g)){
       g<-igraph::as.undirected(g, mode='collapse')
     }

     current_neighbours <- as.integer(igraph::neighbors(g, v=current_node, mode='all'))

     for(neighbour in current_neighbours){
       if(!is.null(prev_node)){
         if(neighbour==prev_node){
           sampling_weights <- c(sampling_weights, igraph::E(g, p=c(current_node, neighbour))$weight / p)
         }else if(igraph::are_adjacent(g, neighbour, prev_node)){
           sampling_weights <- c(sampling_weights, igraph::E(g, p=c(current_node, neighbour))$weight)
         }else{
           sampling_weights <- c(sampling_weights, igraph::E(g, p=c(current_node, neighbour))$weight / q)
         }
       }else{
         sampling_weights <- c(sampling_weights, igraph::E(g, p=c(current_node, neighbour))$weight / q)
       }
     }

     probs <- sampling_weights / sum(sampling_weights)
     next_step <- current_neighbours[sample(length(current_neighbours), size=1, replace=T, prob=probs)]
     return(next_step)
   }

   random_walk <- function(g){
     final_walks <- c()
     vocabulary_size <- length(igraph::V(g))
     nodes <- as.vector(igraph::V(g))


     for(walk_iter in 1:num.walks){
       #message(paste0('Walk number ', walk_iter, ' out of ', num_walks))
       nodes <- sample(nodes)

       #progress_bar <- utils::txtProgressBar(min = 1, max = length(nodes), initial = 1)
       for(node in nodes){
         #utils::setTxtProgressBar(progress_bar, value=1)
         walk <- c(node)

         while(length(walk) < num.steps){
           current_node <- walk[length(walk)]
           if(length(walk) > 1){
             prev_node <- walk[length(walk)-1]
           }else{
             prev_node <- NULL
           }

           next_node <- rw_step(g, prev_node, current_node)
           walk <- c(walk, next_node)
         }
         final_walks[[walk_iter]] <- walk
       }
     }


     return(list(walks = final_walks, vocabulary_size = vocabulary_size))
   }

   create_dataset <- function(walks.list){

     walks <- walks.list$walks
     vocabulary_size <- walks.list$vocabulary_size

     final_targets <- c()
     final_contexts <- c()
     final_labels <- c()

     datasets <- vector(mode = 'list', length = length(walks))

     for(i in 1:length(walks)){
       walk <- walks[[i]]
       skipgrams <- keras::skipgrams(
         walk,
         vocabulary_size = vocabulary_size,
         window_size = window.size,
         negative_samples = num.negative.samples,
         seed = seed
       )

       couples <- skipgrams$couples
       labels <- skipgrams$labels


       for(j in 1:length(couples)){
         pair <- couples[[j]]
         label <- labels[j]
         target <- min(pair)
         context <- max(pair)

         if(target == context){
           next
         }

         final_targets <- c(final_targets, target)
         final_contexts <- c(final_contexts, context)
         final_labels <- c(final_labels, label)
       }

       datasets[[i]] <- data.frame(target = final_targets,
                             context = final_contexts,
                             label = final_labels)
     }


     dataset <- do.call('rbind', datasets)

     return(list(dataset = dataset, vocabulary_size = vocabulary_size))
   }

   node2vec_train <- function(dataset.list){
     dataset <- dataset.list$dataset
     vocabulary_size <- dataset.list$vocabulary_size

     target <- keras::layer_input(name = 'target', shape = c(1))
     context <- keras::layer_input(name = 'context', shape = c(1))

     embed_item <- keras::layer_embedding(input_dim = vocabulary_size + 1,
                                          output_dim = embedding.dim,
                                          embeddings_initializer = 'he_normal',
                                          embeddings_regularizer = keras::regularizer_l2(1e-6),
                                          name = 'embeddings')

     target_embeddings <- embed_item(target)
     context_embeddings <- embed_item(context)
     dot <- keras::layer_dot(axes = 2, normalize = FALSE, name = 'cosine_sim')
     logits <- dot(list(target_embeddings, context_embeddings))
     logits <- keras::k_squeeze(logits, axis = 2)
     model <- keras::keras_model(list(target, context), logits)

     model %>% keras::compile(optimizer = 'adam',
                              loss = keras::loss_binary_crossentropy(from_logits = T),
                              metrics = list('accuracy'))

     embedding_model <- keras::keras_model(inputs = target, outputs = keras::get_layer(model, 'embeddings')$output)

     model %>% keras::fit(x = list(target = dataset$target, context = dataset$context),
                          y = dataset$label,
                          epochs = epochs,
                          batch_size = batch.size)


     return(list(embedding_model = embedding_model, vocabulary_size = vocabulary_size))
   }


   embed_nodes <- function(embedding.model){

     embedding_model <- embedding.model$embedding_model
     vocabulary_size <- embedding.model$vocabulary_size
     out_embeddings <- vector(mode = 'list', length = vocabulary_size)

     for(i in 1:vocabulary_size){
       out_embeddings[[i]] <- stats::predict(embedding_model, i)
     }

     out_embeddings <- do.call('rbind', out_embeddings)

     return(out_embeddings)

   }

   project_embeddings <- function(embeddings){

     if(dim.reduction == 'pca'){
       out <- stats::prcomp(t(embeddings))$rotation[,1:2]

     }else if(dim.reduction == 'tsne'){
       requireNamespace('Rtsne')
       out <- Rtsne::Rtsne(embeddings, perplexity = tsne.perplexity)$Y

     }else if(dim.reduction == 'umap'){
       requireNamespace('umap')
       out <- umap::umap(embeddings)$layout

     }else{
       stop('Dim reduction algorithm is not available!')
     }

     out <- as.data.frame(out)
     colnames(out) <- c('X', 'Y')
     return(out)
   }

   get_graph <- function(tree){
     if(graph.type == 'tree'){
       g <- tree@tree


     }else if(graph.type == 'heterogeneous'){
       g <- tree@heterogeneous

     }else{
       stop('Graph type not recognized/available!')
     }

     #Ensure the weights are inversely proportional to the LV distances (lower weights = further from each other = lower walk step probability)
     igraph::E(g)$weight <- 1/igraph::E(g)$weight

     return(g)
   }

   spectral_embedding <- function(g){

     g <- igraph::as.undirected(g)

     if(embedding.method == 'spectral_adjacency'){
       out <- igraph::embed_adjacency_matrix(g, 2)$X
     }else if(embedding.method == 'spectral_laplacian'){
       out <- igraph::embed_laplacian_matrix(g, 2)$X
     }else{
       stop('Embedding method not recognized!')
     }

     colnames(out) <- c('X', 'Y')

     return(out)
   }

   plot_embeddings <- function(embeddings, tree){

     g <- get_graph(tree)

     sample_id <- paste0(tree@sample_id, collapse = '; ')
     clonotype_id <- tree@clonotype_id

     if(length(clonotype_id) > 5){
       clonotype_id <- 'joined clonotypes'
     }else{
       clonotype_id <- paste0(clonotype_id, collapse = '; ')
     }

     if(is.null(color.by)){
       color.by <- tree@feature_names

       if(is.null(color.by)){
         color.by <- 'node_type'
       }
     }

     final_df <- vector(mode = 'list', length = length(color.by))

     for(i in 1:length(color.by)){
       df <- as.data.frame(embeddings)

       feature <- color.by[i]
       selected_values <- igraph::vertex_attr(g, name = feature)
       if(paste0(feature, '_counts') %in% igraph::vertex_attr_names(g)){
         counts <- igraph::vertex_attr(g, name=paste0(feature, '_counts'))
         max_indices <- lapply(counts, function(x) which.max(x))
         max_values <- unlist(mapply(function(x,y) x[y], selected_values, max_indices))
         feature_values <- max_values
       }else{
         feature_values <- selected_values
       }

       df$feature <- unlist(feature_values)
       df$feature_name <- feature
       df$cell_number <- igraph::vertex_attr(g, name = 'cell_number')
       df$cell_number[is.na(df$cell_number)] <- 1

       final_df[[i]] <- df
     }

     final_df <- do.call('rbind', final_df)

     plot <- ggplot2::ggplot(data = final_df, ggplot2::aes(x = X, y = Y, color = feature)) + #optional: size = cell_number
             ggplot2::geom_point(size = 0.75) +
             ggplot2::theme_bw() +
             ggplot2::theme_classic() +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
             ggplot2::labs(title = paste0(sample_id, ' - ', clonotype_id)) +
             ggplot2::facet_wrap(~feature_name, scales = "free")


     return(plot)
   }

   set.seed(seed)
   
   if(inherits(trees, 'list')){

     plots <- vector(mode = 'list', length = length(trees))
     for(i in 1:length(trees)){
       plots[[i]] <- vector(mode = 'list', length = length(trees[[i]]))

       embeddings <- vector(mode = "list", length = length(trees[[i]]))

       if(parallel){
         requireNamespace('parallel')
         cores <- parallel::detectCores()

         if(embedding.method  == 'node2vec'){
           embeddings <- parallel::mclapply(trees[[i]], mc.cores = cores,
                                                    FUN = function(x) {x %>% get_graph() %>% random_walk() %>% create_dataset() %>% node2vec_train() %>% embed_nodes() %>% project_embeddings()
                                                                       })

         }else if(stringr::str_detect(embedding.method, 'spectral')){
           embeddings <- parallel::mclapply(trees[[i]], mc.cores = cores,
                                                    FUN = function(x) {x %>% get_graph() %>% spectral_embedding()
                                                                       })

         }else{
           stop('Method is not available!')
         }


       }else{
         if(embedding.method == 'node2vec'){
            embeddings <- lapply(trees[[i]], function(x) {x %>% get_graph() %>% random_walk() %>% create_dataset() %>% node2vec_train() %>% embed_nodes() %>% project_embeddings()
                                                               })

         }else if(stringr::str_detect(embedding.method, 'spectral')){
            embeddings <- lapply(trees[[i]], function(x) {x %>% get_graph() %>% spectral_embedding()
                                                         })
         }else{
           stop('Method is not available!')
         }
       }


       for(j in length(trees[[i]])){
         plots[[i]][[j]] <- plot_embeddings(embeddings[[j]], trees[[i]][[j]])
       }

     }


   }else if(inherits(trees, 'AntibodyForests')){

     if(embedding.method == 'node2vec'){
        embeddings <- trees %>% get_graph() %>% random_walk() %>% create_dataset() %>% node2vec_train() %>% embed_nodes() %>% project_embeddings()

     }else if(stringr::str_detect(embedding.method, 'spectral')){
        embeddings <- trees %>% get_graph() %>% spectral_embedding()

     }else{
       stop('Method is not available!')
     }

     plots <- plot_embeddings(embeddings, trees)


   }else{
     stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
   }


   return(plots)
}
