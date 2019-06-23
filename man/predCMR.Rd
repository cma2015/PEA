\name{predCMR}
\alias{predCMR}



\title{
Predicting the CMR status using machine learning-based models
}


\description{
This function is used to perform CMR prediction using machine learning-based
model.
}

\usage{
predCMR(predSeqs, seqLen = 101, motif = "[ag][ag]ac[act]", model, ...)
}


\arguments{
  \item{predSeqs}{a character string specifying the directory of FASTA format sequences.}
  \item{seqLen}{The length of extracted sequence.}
  \item{motif}{A string, which specifies the motif to be searched.}
  \item{model}{A class object representing the model.}
  \item{\dots}{Further parameters passed to predict CMRs. }
}



\value{
A matrix with two columns representing the transcript IDs and corresponding scores.
}


\author{
Chuang Ma, Jingjing Zhai
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
  
  
   predSeq <- system.file("extdata/test_pred.fa", package = "PEA")
   predMat <- predCMR(predSeq, model = psolRes$model)
}
}


\keyword{machine learning}
