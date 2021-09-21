#' @title Simulate B or T cell receptor sequences by variational autoencodes(VAEs) trained with experimental data.
#' @param sequence a vector of seuqnece the model to be trained on
#' @param batch.size set to larger to save time, set to smaller to same computing power
#' @param latent.dim parameter used in VAE model
#' @param intermediate.dim parameter used in VAE model
#' @param epochs parameter used in VAE model
#' @param epsilon.std parameter used in VAE model
#' @param n.train number of sequence to be used in training set, the rest will be in testing set
#' @param n.sample number of new sequence to generate from VAE model
#' @param null.threshold threshold of predicted value to be considered as an existing base, default is 0.05. When generated sequence is too short, lower this threshold.
#' @importFrom dplyr %>%
#' @export

Echidna_vae_generate<-function(sequence,
                               n.train,
                               n.sample,
                               batch.size,
                               latent.dim,
                               intermediate.dim,
                               epochs,
                               epsilon.std,
                               null.threshold){

  if(missing(batch.size)) batch.size<-10L
  if(missing(latent.dim)) latent.dim<-2L
  if(missing(intermediate.dim)) intermediate.dim<-256
  if(missing(epochs)) epochs<-50L
  if(missing(epsilon.std)) epsilon.std<-1
  if(missing(null.threshold)) null.threshold<-0.05
  base<-c("A","T","G","C")

  #require(keras,tensorflow,stats)

  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()


  K <- keras::backend()

  #processing-----------
  sequence<-toupper(sequence)
  data<-.onehot_encoding(sequence,base)
  # Parameters --------------------------------------------------------------

  batch_size <- batch.size
  original_dim <- ncol(data)
  latent_dim <- latent.dim
  intermediate_dim <- intermediate.dim
  epsilon_std <- epsilon.std


  # data into folds
  n_test<-nrow(data)-n.train
  n_fold<-ceiling(nrow(data)/n_test)
  n_overlap<-round((n_fold*n_test-nrow(data))/(n_fold-1))
  cuts<-seq(1,nrow(data),n_test-n_overlap)
  cuts<-c(cuts[1:n_fold],nrow(data))
  all_output<-c()

  for(i in 1:n_fold){
  # Model definition --------------------------------------------------------
  x <- keras::layer_input(shape = c(original_dim))
  h <- keras::layer_dense(x, intermediate_dim, activation = "relu")
  z_mean <- keras::layer_dense(h, latent_dim)
  z_log_var <- keras::layer_dense(h, latent_dim)

  sampling <- function(arg){
    z_mean <- arg[, 1:(latent_dim)]
    z_log_var <- arg[, (latent_dim + 1):(2 * latent_dim)]

    epsilon <- keras::k_random_normal(
      shape = c(keras::k_shape(z_mean)[[1]]),
      mean=0.,
      stddev=epsilon_std
    )

    z_mean + keras::k_exp(z_log_var/2)*epsilon
  }

  # note that "output_shape" isn't necessary with the TensorFlow backend
  z <- keras::layer_concatenate(list(z_mean, z_log_var)) %>%
    keras::layer_lambda(sampling)

  # we instantiate these layers separately so as to reuse them later
  decoder_h <- keras::layer_dense(units = intermediate_dim, activation = "relu")
  decoder_mean <- keras::layer_dense(units = original_dim, activation = "sigmoid")
  h_decoded <- decoder_h(z)
  x_decoded_mean <- decoder_mean(h_decoded)

  # end-to-end autoencoder
  vae <- keras::keras_model(x, x_decoded_mean)

  # encoder, from inputs to latent space
  encoder <- keras::keras_model(x, z_mean)

  # generator, from latent space to reconstructed inputs
  decoder_input <- keras::layer_input(shape = latent_dim)
  h_decoded_2 <- decoder_h(decoder_input)
  x_decoded_mean_2 <- decoder_mean(h_decoded_2)
  generator <- keras::keras_model(decoder_input, x_decoded_mean_2)


  vae_loss <- function(x, x_decoded_mean){
    xent_loss <- (original_dim/1.0)*keras::loss_binary_crossentropy(x, x_decoded_mean)
    kl_loss <- -0.5*keras::k_mean(1 + z_log_var - keras::k_square(z_mean) - keras::k_exp(z_log_var), axis = -1L)
    xent_loss + kl_loss
  }

  vae %>% keras::compile(optimizer = "rmsprop", loss = vae_loss)



  #start rotation-------

  V_range<-c()



    x_test <- data[cuts[i]:cuts[i+1]-1,]
    x_train <- data[setdiff(1:nrow(data),cuts[i]:cuts[i+1]-1),]


    # Model training ----------------------------------------------------------

    vae %>% keras::fit(
      x_train, x_train,
      shuffle = TRUE,
      epochs = epochs,
      batch_size = batch_size,
      validation_data = list(x_test, x_test)
    )


    #Generate output ---------------------------------
    output<-c()
    x_train_encoded<-as.data.frame(stats::predict(encoder, x_train, batch_size = batch_size))

    V_range<-range(c(V_range,x_train_encoded$V1,x_train_encoded$V2))

    sampling.parameter<-list()
    for(s in 1:ceiling(n.sample/n_fold)){
      sampling.parameter[[s]]<-stats::runif(latent.dim,V_range[1],V_range[2])
    }

    for(p in 1:length(sampling.parameter)){
      z_sample<-matrix(sampling.parameter[[p]],ncol = latent_dim)
      matrix<-stats::predict(generator, z_sample)
      matrix<-matrix(matrix,ncol = original_dim/4)
      seq<-.Arg_base(matrix,null.threshold,base)
      output[p]<-seq
    }


    all_output<-c(all_output,output)
  }
  return(all_output)
}


