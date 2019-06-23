\name{runTopGO}
\alias{runTopGO}



\title{
Go enrichment
}


\description{
GO (Gene Ontology) enrichment analysis through "topGO" package.
}


\usage{
runTopGO(geneID, statistic = "fisher", algorithm = "elim",
         topNodes = 20, 
         dataset = "athaliana_eg_gene", plot = TRUE)
}

                  
\arguments{
   \item{geneID}{A character vector representing the gene set to be used for GO enrichment.}
   \item{statistic}{A character string specifing which test to use.}
   \item{algorithm}{A character string specifing which algorithm to use.}
   \item{topNodes}{The number of top GO terms to be included in the results.}
   \item{dataset}{A character string specifing which dataset to use, the default is: athaliana_eg_gene; For other plant species, please use "listDatasets(mart = useMart(biomart = "plants_mart",  host = 'plants.ensembl.org'))" to see all available datasets.}
   \item{plot}{Logical, if TRUE, a dot plot will be generated to represent the significance of each GO term.}
}



\value{
  A list recording enrichment results for three subcategories (MF, BP and CC).
  \item{BP}{The GO enrichment results for BP(biological process).}
  \item{MF}{The GO enrichment results for MF(molecular function).}
  \item{CC}{The GO enrichment results for CC(cellular component).}
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  ##Please ensure packages "biomaRt" and "topGO" are pre-installed.
  library(biomaRt)
  library(topGO)
  mart <- biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", 
                            host = 'plants.ensembl.org')
  GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
  GTOGO <- GTOGO[GTOGO$go_id != '', ]
  all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
  geneID <- sample(all.genes, 100)
  results <- runTopGO(geneID = geneID)
}
}


\keyword{GO enrichment}

