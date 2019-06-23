\name{motifScan}
\alias{motifScan}



\title{
Scanning the motifs in the peak regions
}


\description{
For a given motif, all peak regions will be scanned and then calculates a position weight matrix(PWM) for the analyzed motifs. Finally, a motif logo will be generated to show the positional differences in the frequency of four nucleotides.
}


\usage{
motifScan(sequence, motif = "[ag][ag]ac[act]")
}

                  
\arguments{
   \item{sequence}{A character vector representing the file name of the sequence in FASTA format.}
   \item{motif}{A string, which specifies the motif to be searched.}
}



\value{
  A PWM matrix and motif logo.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  testSeq <- system.file("extdata/test.fa", package = "PEA")
  results <- motifScan(sequence = testSeq)
}
}


\keyword{CMR Annotation}

