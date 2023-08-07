#' This PlatypusML_classification function takes as input encoded features obtained using the PlatypusML_extract_features function. The function runs cross validation on a specified number of folds
#' for different classification models and reports the AUC scores and ROC curves.
#' @param matrix Matrix. Output of the PlatypusML_extract_features function, with the last column storing the label.
#' @param label.1 String. The label of the first class.
#' @param label.2 String. The label of the second class.
#' @param proportion. Vector of size 2 (floats between 0 and 1 that need to sum up to 1). Specifies the proportions for the two classes.
#' The smaller proportion will be assigned to the minority class by default.
#' @param random.seed Integer. The seed to be set when sampling for balancing the dataset.
#' @return This function returns a matrix containing equal number of samples for the two classes.
#' @examples
#' \dontrun{
#' PlatypusML_balance(matrix, label.1 = a, label.2 = b)
#' }




PlatypusML_balance <- function(matrix,
                               label.1,
                               label.2,
                               proportion,
                               random.seed){


  if(missing(matrix)) stop("Please provide matrix input for this function")
  if(missing(label.1)) stop("Please provide the label of the first class")
  if(missing(label.2)) stop("Please provide the label of the second class")
  if(missing(proportion)) proportion <- c(0.5, 0.5)
  if(missing(random.seed)) random.seed = 1

  set.seed(random.seed)
  class1 <- subset(matrix, matrix[,ncol(matrix)]==label.1)
  class2 <- subset(matrix, matrix[,ncol(matrix)]==label.2)

  #assign the different proportions
  prop1 <- min (proportion)
  prop2 <- max (proportion)

  #if class 1 is the minority
  if (nrow(class1)<= nrow(class2)){


    size2 <- floor(prop2/prop1*nrow(class1))

    #if the proportions are not suitable, throw error
    if (size2>nrow(class2)){
      stop("Minimum proportion is too low, please increase it.")
    }
    class2 <- class2[sample.int(nrow(class2), size2),]
  }  else {

    size1 <- floor(prop1/prop2*nrow(class2))

    if (size1>nrow(class1)){
      stop("Minimum proportion is too low, please increase it.")
    }
    class1 <- class1[sample.int(nrow(class1), size1),]
  }


  balanced_matrix <- rbind(class1, class2)
  shuffled_balanced_matrix <-  balanced_matrix[sample(1:nrow( balanced_matrix)),]

  return (shuffled_balanced_matrix)

}
