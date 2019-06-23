\name{extractSeqs}
\alias{extractSeqs}



\title{
Extract sequences according to samples
}


\description{
Extracting the sequences of fixed length according to the samples.
}


\usage{
extractSeqs(RNAseq, samples, seqLen = 101)
}

                  
\arguments{
   \item{RNAseq}{A character string specified the file name of cDNA sequences.}
   \item{samples}{A character vetor specified the samples.}
   \item{seqLen}{The length of extracted sequence.}
}



\value{
  A list containing the extracted sequences.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{

  ##Load data
  data(sampleData)
  positiveSamples <- sampleData$positives
  cDNA <- system.file("extdata/chr1_cdna.fa", package = "PEA")
  posSeq <- extractSeqs(RNAseq = cDNA, samples = positiveSamples, seqLen = 101)
}
}



