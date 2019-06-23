\name{motifDetect}
\alias{motifDetect}



\title{
De novo motif discovery in the peak regions
}


\description{
The de novo motif discovery is performed through genetic algorithm in 'rGADEM' package.
}


\usage{
motifDetect(sequence, plot = T, ...)
}

                  
\arguments{
   \item{sequence}{A character vector representing the file name of the sequence in FASTA format.}
   \item{plot}{Logical, where TRUE indicates that perform visualization for the results.}
   \item{\dots}{Further paramters used for motif detection, see GADEM.}
}



\value{
  A list containing consensus motifs in PWM (position weight matrix) format.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  testSeq <- system.file("extdata/test.fa", package = "PEA")
  results <- motifDetect(sequence = testSeq)
}
}


\keyword{motif detection}

