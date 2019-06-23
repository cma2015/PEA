\name{findConfidentPosSamples}
\alias{findConfidentPosSamples}



\title{
Find confident positive samples
}


\description{
The confident positive samples are identified if only one CMR signal existed in the corresponding peak.
}


\usage{
findConfidentPosSamples(peaks, RNAseq = NULL, motifPos = NULL, ...)
}

                  
\arguments{
   \item{peaks}{A BED format matrix containing the peak regions.}
   \item{RNAseq}{A character string specified the file name of cDNA sequences.}
   \item{motifPos}{A list containing the position of the consensus motif for each transcript.}
   \item{\dots}{Further paramters used for find confident positive samples.}
}



\value{
  A list containing two elements: the positive samples and its transcript IDs.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  
  ##Load data
  input.bam <- system.file("extdata/chr1_input_test.bam", package = "PEA")
  RIP.bam <- system.file("extdata/chr1_RIP_test.bam", package = "PEA")
  refGenome <- system.file("extdata/chromosome1.fa", package = "PEA")
  cDNA <- system.file("extdata/chr1_cdna.fa", package = "PEA")
  GTF <- system.file("extdata/chromosome1.gtf", package = "PEA")

  ##m6A peak calling using sliding window-based method
  cmrMat <- CMRCalling(CMR = "m6A", IPBAM = RIP.bam, inputBAM = input.bam,
                       method = "SlidingWindow", mappedInput = 17472, 
                       mappedRIP = 20072, refGenome = refGenome)


  #Convert genomic position to cDNA position
  peaks <- G2T(bedPos = cmrMat, GTF = GTF)
  
  ##Search consensus motif in cDNA sequence
  motifPos <- searchMotifPos(sequence = cDNA)
  
  ##Find confident positive samples
  posSamples <- findConfidentPosSamples(peaks = peaks, motifPos = motifPos)
  
}
}



