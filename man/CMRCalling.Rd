\name{CMRCalling}
\alias{CMRCalling}



\title{
  Identify CMRs (chemical modifications of RNA) from epitranscriptome sequencing data
}


\description{
  CMRCalling is designed to identify different types of CMRs from epitranscriptome sequencing data.
}


\usage{
  CMRCalling(CMR = c("m6A", "m6Am", "m5C", "hm5C", "pseudoU", "m1A"),
             cpus = 1, IPBAM = NULL, inputBAM = NULL, 
             GTF = NULL,  paired = FALSE, ...)
}

                  
\arguments{
  \item{CMR}{A string, which specifies the CMR type.}
  \item{cpus}{A integer number specifying the number of cpus to be used for parallel computing.}
  \item{IPBAM}{A string vector of file name, which specifies the IP samples in BAM format.}
  \item{inputBAM}{A string vector of file name, which specifies the input control samples in BAM format.}
  \item{GTF}{A string vector of file name, which specifies the reference genome annotation in GTF format.}
  \item{paired}{Logical, if TRUE, the paired-end parameters will be used for CMR calling, otherwise, the single-end parameters would be used.}
  \item{\dots}{Further parameters used CMR calling.}
}



\value{
  A matrix in BED format containing the CMRs.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  
  #Preparing sample data for peak calling
  input.bam <- system.file("extdata/chr1_input_test.bam", package = "PEA")
  RIP.bam <- system.file("extdata/chr1_RIP_test.bam", package = "PEA")
  refGenome <- system.file("extdata/chromosome1.fa", package = "PEA")
  GTF <- system.file("extdata/chromosome1.gtf", package = "PEA")
  
  
  #m6A peak calling using sliding window-based method
  #mappedInput: represents the number of mapped reads number in the input control samples.
  #mappedRIP: represents the number of mapped reads number in the RIP samples.
  cmrMat <- CMRCalling(CMR = "m6A", IPBAM = RIP.bam, inputBAM = input.bam,
                       method = "SlidingWindow", mappedInput = 17472, 
                       mappedRIP = 20072, refGenome = refGenome)
                       
  #m6A peak calling using exomePeak
  results <- CMRCalling(CMR = "m6A", method = "exomePeak", IPBAM = RIP.bam,
                       inputBAM = input.bam, GTF = GTF)
                      
  #m6A peak calling using BayesPeak
  results <- CMRCalling(CMR = "m6A", method = "BayesPeak", IPBAM = RIP.bam,
                       inputBAM = input.bam, GTF = GTF)  
 
  #m6A peak calling using MetPeak
  results <- CMRCalling(CMR = "m6A", method = "MetPeak", IPBAM = RIP.bam,
                       inputBAM = input.bam, GTF = GTF)
                       
  #m6A peak calling using MACS2
  results <- CMRCalling(CMR = "m6A", method = "MACS2", IPBAM = RIP.bam,
                       inputBAM = input.bam, GTF = GTF, ...="--nomodel")
                       
  #Pseudo uridine ratio
  inputBAM <- system.file("extdata/pseudoU_test.bam", package = "PEA")
  refGenome <- system.file("extdata/chr1_cdna.fa", package = "PEA")
  results <- CMRCalling(CMR = "pseudoU", inputBAM = inputBAM,
                        cpus = 4, refGenome = refGenome)
                          
  #m5C calling
  fq <- system.file("extdata/test.fq", package = "PEA")
  refGenome <- system.file("extdata/chromosome1.fa", package = "PEA")
  results <- CMRCalling(CMR = "m5C", fq = fq, GTF = GTF,
                        cpus = 4, refGenome = refGenome)
                          
}
}




\keyword{CMR calling, epitranscriptome}

