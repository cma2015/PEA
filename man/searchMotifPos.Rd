\name{searchMotifPos}
\alias{searchMotifPos}



\title{
Scanning the motif in the sequence
}


\description{
For a given motif, all sequence will be scanned return the central position for the motif.
}


\usage{
searchMotifPos(sequence, motif = "[ag][ag]ac[act]", cenPos = 2)
}

                  
\arguments{
   \item{sequence}{A character vector representing the file name of the sequence in FASTA format.}
   \item{motif}{A string, which specifies the motif to be searched.}
   \item{cenPos}{Integer represents the central position for the motif.}
}



\value{
  A list containing the motif central position.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{

  cDNA <- system.file("extdata/chr1_cdna.fa", package = "PEA")
  results <- searchMotifPos(sequence = cDNA)
}
}


\keyword{motif}

