\name{QCBAM}
\alias{QCBAM}



\title{
  Quality control for BAM format files
}


\description{
  Performing quality control for mapping results in BAM format.
}


\usage{
  QCBAM(BAM, quality = c(35, 30), paired = F)
}

                  
\arguments{
  \item{BAM}{A character vector representing file name in bam format.}
  \item{quality}{A vector which specifies the base-call quality score to be filtered for single-end and paired-end, respectively.}
  \item{paired}{Logical, where TRUE indicates the BAM files are generated with paired-end reads.}
}



\value{
  A list containing filtered files which passed the quality control.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
   BAM <- system.file("extdata/chr1_input_test.bam", package = "PEA")
   results <- QCBAM(BAM = BAM)
}
}




