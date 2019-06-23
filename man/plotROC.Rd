\name{plotROC}
\alias{plotROC}



\title{
Plotting receiver operating characteristic(ROC) Curves
}


\description{
This function plots ROC curves for estimating the performance of machine learning-based 
classification model in cross validation experiments.
}


\usage{
plotROC(cvRes)
}

\arguments{
  \item{cvRes}{results from the "cross_validation" function.}
}



\value{
 A ROC plot
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
   
   ##Perform five-fold cross-validation using random forest algorithm
   cvRes <- cross_validation(method = "randomForest", featureMat = featureMat,
                             positives = positiveSamples, negatives = negativeSamples, 
                             cross = 5, cpus = 1)
  plotROC(cvRes = cvRes)
                          
}
}


\keyword{ machine learning }

