\name{pseudoURatio}
\alias{pseudoURatio}



\title{
Pseudouridine ratio
}


\description{
  Calculate the pseudouridine ratio from the pseudo-seq data.
}


\usage{
  pseudoURatio(refGenome, inputBAM, cpus = 1)
}

                  
\arguments{
   \item{refGenome}{A character vector representing file name of reference genome in FASTA format.}
  \item{inputBAM}{A string vector of file name, which specifies the input control samples in BAM format.}
   \item{cpus}{an integer number specifying the number of cpus to be used for parallel computing.}
}




\value{
  A list  containing the position and ratio for each pseudouridine.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{

  inputBAM <- system.file("extdata/pseudoU_test.bam", package = "PEA")
  referenceGenome <- system.file("extdata/chr1_cdna.fa", package = "PEA")
  results <- pseudoURatio(refGenome = referenceGenome,
                          inputBAM = inputBAM, cpus = 4)
}
}

\references{Carlile, T.M., et al. Pseudouridine profiling reveals regulated mRNA pseudouridylation in yeast and human cells. Nature 2014;515(7525):143-146.}


