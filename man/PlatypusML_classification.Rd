% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlatypusML_classification.R
\name{PlatypusML_classification}
\alias{PlatypusML_classification}
\title{This PlatypusML_classification function takes as input encoded features obtained using the PlatypusML_extract_features function. The function runs cross validation on a specified number of folds
for different classification models and reports the AUC scores and ROC curves.}
\usage{
PlatypusML_classification(features, cv.folds, balancing, proportion)
}
\arguments{
\item{features}{Matrix. Output of the PlatypusML_extract_features function, containing the desired label in the last column.}

\item{cv.folds}{Integer. The number of folds to be used in cross validation}

\item{balancing}{Boolean. Whether to perform class balancing. Defaults to TRUE.}

\item{proportion.}{Vector of size 2 (floats between 0 and 1 that need to sum up to 1). Specifies the proportions for the two classes.
The smaller proportion will be assigned to the minority class by default. Defaults to c(0.5,0.5).}
}
\value{
This function returns a list containing [["combined"]] summary plot with ROC & confusion matrices, [["ROC"]] the ROC curve, [["confusion"]] confusion matrices for each classifier.
}
\description{
This PlatypusML_classification function takes as input encoded features obtained using the PlatypusML_extract_features function. The function runs cross validation on a specified number of folds
for different classification models and reports the AUC scores and ROC curves.
}
\examples{
\dontrun{
To classify and obtain the performance of different models, using extracted and encoded features.

#extract features
features_VDJ_GP33_binder <- PlatypusML_feature_extraction_VDJ(VGM = VGM,
which.features = c("VDJ_cdr3s_nt"),
which.encoding = c("kmer"),
parameters.encoding.nt = c(3),
which.label = "GP33_binder")

#classify
classifier_GP33_binder <- classification(features = features_VDJ_GP33_binder,
cv.folds = 5,
balancing = TRUE)

#view summary
classifier_GP33_binder$combined
}
}
