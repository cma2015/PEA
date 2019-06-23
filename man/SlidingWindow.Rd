\name{SlidingWindow}
\alias{SlidingWindow}



\title{
Fisher exact test based sliding window
}


\description{
  Peak calling through Fisher exact test based sliding window method (see References).
}


\usage{
  SlidingWindow(input, RIP, mappedInput = NULL, 
                mappedRIP = NULL, level = 0.05, 
                ratio = -1, cpus = 1, readsCount = 10,
                concatenate = 4, ...)
}

                  
\arguments{
  \item{input}{A string vector of file name, which specifies the input control samples in BAM format.}
    \item{RIP}{A string vector of file name, which specifies the IP samples in BAM format.}
    \item{mappedInput}{Integer, which represents the number of mapped reads number in the input control samples.}
    \item{mappedRIP}{Integer, which represents the number of mapped reads number in the IP samples.}
    \item{level}{Significance level for the p-value.}
    \item{ratio}{The ratio between the normalized mapped reads number in IP and input samples, respectively.}
    \item{cpus}{The number of threads used for parallel computing.}
    \item{readsCount}{An interger specified the minimum number of reads count in IP sample.}
    \item{concatenate}{The minimum number of continuous window to be concatenated.}
    \item{\dots}{Further parameters used for peak calling.}
}




\value{
  A matrix in BED format containing the enriched peaks.
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
  
  
  ########################peak calling##############################
  peaks <- SlidingWindow(RIP = RIP.bam, input = input.bam,  
                         mappedInput = 17472, 
                         mappedRIP = 20072, refGenome = refGenome)
                          
}
}




\keyword{peak calling, epitranscriptome}

\references{Meyer, K.D., ea al. Comprehensive analysis of mRNA methylation reveals enrichment in 3'UTRs and near stop codons. Cell 2012; 149(7):1635-1646}
