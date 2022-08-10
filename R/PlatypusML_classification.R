#' This PlatypusML_classification function takes as input encoded features obtained using the PlatypusML_extract_features function. The function runs cross validation on a specified number of folds 
#' for different classification models and reports the AUC scores and ROC curves. 
#' @param features Matrix. Output of the PlatypusML_extract_features function, containing the desired label in the last column.
#' @param cv.folds Integer. The number of folds to be used in cross validation
#' @param balancing Boolean. Whether to perform class balancing. Defaults to TRUE. 
#' @param proportion. Vector of size 2 (floats between 0 and 1 that need to sum up to 1). Specifies the proportions for the two classes.
#' The smaller proportion will be assigned to the minority class by default. Defaults to c(0.5,0.5).
#' @return This function returns a list containing [["combined"]] summary plot with ROC & confusion matrices, [["ROC"]] the ROC curve, [["confusion"]] confusion matrices for each classifier.
#' #' @examples
#' \dontrun{
#' To classify and obtain the performance of different models, using extracted and encoded features.
#' 
#' #extract features
#' features_VDJ_GP33_binder <- PlatypusML_feature_extraction_VDJ(VGM = VGM,
#' which.features = c("VDJ_cdr3s_nt"), 
#' which.encoding = c("kmer"),
#' parameters.encoding.nt = c(3),
#' which.label = "GP33_binder")
#' 
#' #classify
#' classifier_GP33_binder <- classification(features = features_VDJ_GP33_binder,
#' cv.folds = 5,
#' balancing = TRUE)
#' 
#' #view summary
#' classifier_GP33_binder$combined
#' }


PlatypusML_classification <- function(features,
                     cv.folds,
                     balancing, 
                     proportion){
  
  
  ###### include functions and libraries ############
  source("PlatypusML_balance.R")
  
  ###### setting defaults ########
  if(missing(cv.folds)) cv.folds <-5
  if(missing(balancing)) balancing <-TRUE
  if(missing(proportion)) proportion <- c(0.5,0.5)
  #balance the data if requested
  if(balancing){
      label.1 <- unique(features[,ncol(features)])[1]
      label.2 <- unique(features[,ncol(features)])[2]
      features<- PlatypusML_balance(features, label.1 =  label.1, label.2 = label.2, proportion = proportion)
    }

  ggplot2::theme_set(ggplot2::theme_bw())
  
  features[,ncol(features)] <- as.numeric(as.factor(features[,ncol(features)])) - 1
  
  #creating folds 
  CV_folds <- caret::createFolds(features[,ncol(features)], k = cv.folds)
  probs<- NULL
  y_tests <- NULL
  for (i in CV_folds){
    
    models <- list()
    
    x_train <- as.matrix(features[-i,-ncol(features)])
    x_test <- as.matrix(features[i, -ncol(features)])
    y_train <-  unlist(features[-i,ncol(features)])
    y_test <-  unlist(features[i,ncol(features)])
    
    models[['svm']] <- e1071::svm(formula = y_train ~ .,
                                  data = x_train,
                                  kernel = 'radial',
                                  probability = TRUE)
    models[['xgb']] <- xgboost::xgb.train(params = list(
      booster = "gbtree",
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = 8),
      data = xgboost::xgb.DMatrix(
        label = y_train,
        data = x_train),
      nrounds = 25)
    models[["gnb"]] <- naivebayes::naive_bayes(
      x = x_train,
      y = as.factor(y_train))
    
    models[["rForest"]] <- randomForest::randomForest(
      x = x_train,
      y = as.factor(y_train))
    
    logreg_formula <- paste(colnames(features)[ncol(features)], " ~ .", sep = "")
    models[['logreg']] <- glm(formula = logreg_formula,
                              data = features[-i,],
                              family = "binomial")
    
    #extract the current probabilities and the y tests for each of the folds 
    current_probs <- data.frame(sapply(models, function (xx) {
      if ((substr(xx$call, 1, 4) =="glm") [1]){
        x_test_ap <- features[i,]
        return(predict(xx, newdata = x_test_ap))
      }
      else{
        x_test_ap <- x_test
      }
      
      pred <- predict(xx, newdata = x_test_ap, type = "prob", probability = TRUE)
      if (ncol(data.frame(pred))>1){
        pred <- data.frame(pred[,2])
      }
      if ((substr(xx$call, 1, 4) =="SVM") [1]){
        pred <- attr(pred, "probabilities")
      }
      return (pred)
    }))
    
    colnames(current_probs) <- names(models)
    rownames(current_probs) <- i
    probs<- rbind(probs, current_probs)
    y_tests<- c(y_tests, y_test)
  }
  
  # ROC curve
  rocs <- lapply(probs, function(x) {
    pROC::roc(y_tests, x, ci = T)
  })
  
  # Plotting ROC
  ROC.plot <- pROC::ggroc(rocs) +
    ggplot2::geom_abline(slope=1, linetype = "dashed", color = "grey", intercept = 1) +
    ggplot2::ggtitle(paste0("ROC curves: ", colnames(features[ncol(features)]))) + 
    ggplot2::labs(color = "model") +
    ggplot2::theme_classic() +
    ggplot2::annotate(geom = "text", x = 0.25, y = seq(0.1, 0.4, length.out = length(rocs)),size=3, 
             label = paste("auc", names(rocs), ":", sapply(rocs, function(x) as.character(round(x$auc, 3)))))
  
  # confusion matrix for each model
  #convert probabilities into binary labels 0s and 1s
  probs$svm[probs$svm>=0.5] <- 1
  probs$svm[probs$svm<0.5] <- 0
  probs$svm <- factor(probs$svm, levels=c("0","1"))
  probs$logreg[probs$logreg>=0.0] <- 1
  probs$logreg[probs$logreg<0.0] <- 0
  probs$logreg <- factor(probs$logreg, levels=c("0","1"))
  probs$xgb[probs$xgb>=0.5] <- 1
  probs$xgb[probs$xgb<0.5] <- 0
  probs$xgb <- factor(probs$xgb, levels=c("0","1"))
  probs$gnb[probs$gnb>=0.5] <- 1
  probs$gnb[probs$gnb<0.5] <- 0
  probs$gnb <- factor(probs$gnb, levels=c("0","1"))
  probs$rForest [probs$rForest>=0.5] <- 1
  probs$rForest [probs$rForest<0.5] <- 0
  probs$rForest <- factor(probs$rForest, levels=c("0","1"))
  truth_estimate <- probs
  truth_estimate$truth <- factor(y_tests, levels=c("0","1"))

  confusion.matrix <- list()
  model_names <- colnames(truth_estimate)
  model_names <- model_names[-length(model_names)]
  for (i in model_names){
    cm <- yardstick::conf_mat(truth_estimate, truth=truth, estimate=i)
    #from https://stackoverflow.com/questions/37897252/plot-confusion-matrix-in-r-using-ggplot
    table <- data.frame(cm$table)
    
    plotTable <- table %>%
      dplyr::mutate(classification = ifelse(table$Prediction == table$Truth, "correct", "incorrect"))
    
    conf.plot <- ggplot2::ggplot(data = plotTable, mapping = ggplot2::aes(x = Truth, y = Prediction, fill = classification, alpha = Freq)) +
      ggplot2::geom_tile(show.legend = FALSE) +
      ggplot2::geom_text(ggplot2::aes(label = Freq), vjust = .5, alpha = 1) +
      ggplot2::scale_fill_manual(values = c(correct = "green", incorrect = "red")) +
      ggplot2::theme_classic() + ggplot2::ggtitle(i) +
      ggplot2::xlim(rev(levels(table$Truth))) + 
      ggplot2::labs(x="true", y="predicted")
    
    confusion.matrix <- append(confusion.matrix, list(conf.plot))
  }
  
  #combining ROC plot & the 5 confusion matrices
  final.plot <- ggpubr::ggarrange(
    ROC.plot, 
    ggpubr::ggarrange(confusion.matrix[[1]], confusion.matrix[[2]], confusion.matrix[[3]], confusion.matrix[[4]], confusion.matrix[[5]], nrow=2, ncol=3, labels=c("B", "C", "D", "E", "F")),
    ncol=2,
    labels="A"
  )
  
  return(list(ROC=ROC.plot, confusion=confusion.matrix, combined = final.plot))
}