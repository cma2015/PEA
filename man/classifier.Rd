\name{classifier}
\alias{classifier}



\title{
Building machine learning-based classification models
}


\description{
This function builds classification models with different machine learning algorithms including
 random forest (randomForest) and support vector machine (svm)
}

\usage{
classifier( method = c("randomForest", "svm"), featureMat,
            positiveSamples, negativeSamples, ...)
}



\arguments{
  \item{method}{a character string specifying machine learning method used to construct
  prediction models. Currently the user-optional values are "randomForest" and "svm".}
  \item{featureMat}{a numeric matrix recording feature values for each sample.}
  \item{positiveSamples}{a character vector recording positive samples used to
   train classification model.}
  \item{negativeSamples}{a character vector recording negative samples used to
   train classification model. }
  \item{\dots}{Further parameters passed to train classification model. }
}


\details{

  For the random forest algorithm, the important parameters are mtry (number of features randomly
   selected for bulding the decision tree. Default:sqrt(ncol(featureMat))) and ntree (number of 
   trees to be built. Default:500). More information about the parameters related to random forest
    can be found in R package RandomForest.

  For the svm algorithm, the default kernel function is the "radial"(radial basis function; RBF).
   Other optional kernels are the "linear", the "polynomial", and the "sigmoid" fuctions. The range
    of degree (for the kernel type of polynomial), gamma (for radial), coef0 (for polynomial 
    and sigmoid) and cost (the cost parameters for all kernels) should be designated. Please refer
   the description of "svm" function in R package e1071 for more information 
   about the parameters of svm. 
}


\value{
An object of class model.
}


\author{
Chuang Ma, Jingjing Zhai
}

\examples{

\dontrun{

   ##Load samples and feature matrix
   data(sampleData)
   
   positiveSamples <- sampleData$positives
   negativeSamples <- sampleData$negatives
   featureMat <- sampleData$featureMat
   
   ##Using ramdom forest algorithm to build machine learning model
   cl <- classifier( method = "randomForest", featureMat = featureMat, 
                     positiveSamples = positiveSamples, 
                     negativeSamples = negativeSamples)
}
}


\keyword{machine learning}
