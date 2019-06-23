\docType{data}
\name{sampleData}
\alias{sample data}



\title{
 example data for PEA
}


\description{
positive and negative samples (m6A/non-m6A modifications) and sequence-based
feature matrix.
}

\usage{
data(sampleData)
}



\examples{
\dontrun{
   ##load datasets
   data(sampleData)
   
   ##load positive samples
   positives <- PEA$positives
   
   ##load negative samples
   negatives <- PEA$negatives
   
   ##Load sequence-based feature matrix
   featuresMat <- PEA$featureMat
 }
}


\keyword{machine learning}
