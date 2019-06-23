\name{PSOL}
\alias{PSOL}


\title{
Machine learning-based PSOL(positive samples only learning algorithm) algorithm
}


\description{
Firstly, this function selects an initial negative set with the machine learning(ML)-based positive-only sample learning (PSOL) algorithm. The PSOL algorithm has been previously applied to predict stress-response genes and genomic loci encoding functional noncoding RNAs (see References).
}


\usage{
PSOL(featureMatrix = NULL, positives, balanced = FALSE, ratio = 10, 
     unlabels, cpus = 1, PSOLResDic = NULL, iterator = 10, TPR = 0.995,
     cross = 5, method = c("randomForest", "svm"), ... )
}




\arguments{
  \item{featureMatrix}{A numeric matrix recording the features for all sample.}
  \item{positives}{A character vector recording positive samples}
  \item{balanced}{A logical value, where TRUE represents balance the positive and negative samples according to the ratio.}
  \item{ratio}{A numuric value of the the ratio between negative and positive samples.}
  \item{unlabels}{A character vector recording unlabeled samples.}
  \item{cpus}{An integer number specifying the number of cpus will be used for parallel computing.}
  \item{PSOLResDic}{A character string, PSOL Result directory}
  \item{iterator}{An integer calue, iterator times.}
  \item{TPR}{A numeric value ranged from 0 to 1.0, used to select the prediction score cutoff.}
  \item{cross}{An integer value, cross-times cross validation.}
  \item{method}{A character string, machine learing method.}
  \item{\dots}{Further parameters used for PSOL.}
}



\value{
A list containing three components:
  \item{positives }{a character vector including the input positive samples.}
  \item{negatives }{a character vector recording the selected negative samples.}
  \item{model }{a class object representing the model.}
}



\author{
Jingjing Zhai, Chuang Ma, Jie Song.
}

\references{
[1] Chunlin Wang, Chris Ding, Richard F. Meraz and Stephen R. Holbrook. PSoL: a positive sample 
only learning algorithm for finding non-coding RNA genes. Bioinformatics, 2006, 22(21): 2590-2596.

[2] Chuang Ma, Xiangfeng Wang. Machine learning-based differential network analysis: a case study
 of stress-responsive transcriptomes in Arabidopsis thaliana. 2013(Submitted).
}


\examples{

\dontrun{
     
   data(sampleData)
   
   positiveSamples <- sampleData$positives
   negativeSamples <- sampleData$negatives
   featureMat <- sampleData$featureMat
   
   #Saving psol results in working directory
   dir.create("psol", showWarnings = F)
   psolResDic <- paste0(getwd(), "/psol/")
   psolRes <- PSOL(featureMatrix = featureMat, positives = positiveSamples,
                   unlabels = negativeSamples, PSOLResDic = psolResDic,
                   cpus = 1)

}

}



\keyword{ machine learning }
\keyword{ PSOL }
