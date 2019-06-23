\name{G2T}
\alias{G2T}



\title{
  Convert the genomic position to transcriptomic position
}


\description{
According to the reference genome GTF annotation file, transform the genomic position to transcriptomic position, if one position can be mapped multiple transcripts, all of them will be retained.
}


\usage{
G2T(bedPos, GTF)
}

                  
\arguments{
  \item{bedPos}{A matrix in BED format containing genomic position.}
   \item{GTF}{A string vector of file name, which specifies the reference genome annotation in GTF format.}
}



\value{
  A list recording transcriptomic position.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  
 #######Load data
 input.bam <- system.file("extdata/chr1_input_test.bam", package = "PEA")
 RIP.bam <- system.file("extdata/chr1_RIP_test.bam", package = "PEA")
 refGenome <- system.file("extdata/chromosome1.fa", package = "PEA")
 GTF <- system.file("extdata/chromosome1.gtf", package = "PEA")


 cmrMat <- CMRCalling(CMR = "m6A", IPBAM = RIP.bam, inputBAM = input.bam,
                      method = "SlidingWindow", mappedInput = 17472, 
                      mappedRIP = 20072, refGenome = refGenome)

 ###Convert genomic position to cDNA position
 peaks <- G2T(bedPos = cmrMat, GTF = GTF)
                          
}
}




\keyword{CMR calling}

