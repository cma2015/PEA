

#' @export
runTopGO <- function(geneID, statistic = "fisher", algorithm = "elim",
                     topNodes = 20, 
                     dataset = "athaliana_eg_gene", plot = TRUE){
  
  # if(!require(biomaRt)){
  #   source("https://bioconductor.org/biocLite.R")
  #   biocLite("biomaRt")
  # }
  # 
  # if(!require(topGO)){
  #   source("https://bioconductor.org/biocLite.R")
  #   biocLite("topGO")
  # }
  
  mart <- useMart(biomart = "plants_mart", dataset = dataset, host = 'plants.ensembl.org')
  GTOGO <- getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
  GTOGO <- GTOGO[GTOGO$go_id != '', ]
  geneID2GO <- by(GTOGO$go_id, GTOGO$ensembl_gene_id, function(x) as.character(x))
  
  all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
  int.genes <- geneID
  int.genes <- intersect(int.genes, names(geneID2GO))
  int.genes <- factor(as.integer(all.genes %in% int.genes))
  names(int.genes) = all.genes
  
  go.obj.BP <- new("topGOdata", ontology='BP',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GO)
  
  go.obj.MF <- new("topGOdata", ontology='MF',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GO)
  
  go.obj.CC <- new("topGOdata", ontology='CC',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GO)
  
  ##########retrieve the gene list related to a GO ID######################
  allGO.BP <- genesInTerm(object = go.obj.BP)
  allGO.MF <- genesInTerm(object = go.obj.MF)
  allGO.CC <- genesInTerm(object = go.obj.CC)
  
  #########retrive the significant GO terms
  results.BP <- runTest(go.obj.BP, algorithm = algorithm, statistic = statistic)
  results.tab.BP <- GenTable(object = go.obj.BP, elimFisher = results.BP,
                             topNodes = topNodes)
  gene.BP <- genesInTerm(object = go.obj.BP, whichGO = results.tab.BP$GO.ID)
  inter.gene.BP <- lapply(X = gene.BP, FUN = function(x) intersect(x, geneID))
  inter.gene.BP <- unlist(lapply(inter.gene.BP, function(x) paste(x, collapse = ';')))
  results.tab.BP$significantGene <- inter.gene.BP
  
  if(length(which(results.tab.BP$elimFisher == "< 1e-30")) != 0){
    results.tab.BP[which(results.tab.BP$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  }
  
  
  results.MF <- runTest(go.obj.MF, algorithm = algorithm, statistic = statistic)
  results.tab.MF <- GenTable(object = go.obj.MF, elimFisher = results.MF, 
                             topNodes = topNodes)
  gene.MF <- genesInTerm(object = go.obj.MF, whichGO = results.tab.MF$GO.ID)
  inter.gene.MF <- lapply(X = gene.MF, FUN = function(x) intersect(x, geneID))
  inter.gene.MF <- unlist(lapply(inter.gene.MF, function(x) paste(x, collapse = ';')))
  results.tab.MF$significantGene <- inter.gene.MF
  if(length(which(results.tab.MF$elimFisher == "< 1e-30")) != 0){
    results.tab.MF[which(results.tab.MF$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  }
  
  results.CC <- runTest(go.obj.CC, algorithm = algorithm, statistic = statistic)
  results.tab.CC <- GenTable(object = go.obj.CC, elimFisher = results.CC, 
                             topNodes = topNodes)
  gene.CC <- genesInTerm(object = go.obj.CC, whichGO = results.tab.CC$GO.ID)
  inter.gene.CC <- lapply(X = gene.CC, FUN = function(x) intersect(x, geneID))
  inter.gene.CC <- unlist(lapply(inter.gene.CC, function(x) paste(x, collapse = ';')))
  
  results.tab.CC$significantGene <- inter.gene.CC
  if(length(which(results.tab.CC$elimFisher == "< 1e-30")) != 0){
    results.tab.CC[which(results.tab.CC$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  }
  
  
  if(plot){
    df <- data.frame(Category = c(rep("BP", topNodes), rep("CC", topNodes), rep("MF", topNodes)), 
                     x = c(results.tab.BP$Significant, results.tab.CC$Significant, 
                           results.tab.MF$Significant),
                     y = c(-log10(as.numeric(results.tab.BP$elimFisher)), 
                           -log10(as.numeric(results.tab.CC$elimFisher)), 
                           -log10(as.numeric(results.tab.MF$elimFisher))),
                     size = c(-log10(as.numeric(results.tab.BP$elimFisher)),
                              -log10(as.numeric(results.tab.CC$elimFisher)), 
                              -log10(as.numeric(results.tab.MF$elimFisher)))
    )
    
    kk <- ggplot(data = df, aes(x = x, y = y)) + 
      geom_point(aes(color = Category, size = size)) + 
      scale_size_continuous(range = c(2,10)) + 
      labs(x = "The number of significant genes", y = "The adjusted p-values for each GO term")
    print(kk)
  }
  
  results <- list(BP = results.tab.BP, CC = results.tab.CC, MF = results.tab.MF)
  results
}



#' @export
searchMotifPos <- function(sequence, motif = "[ag][ag]ac[act]", cenPos = 2){
  Seqs <- read.fasta(sequence, as.string = T)
  if(length(Seqs) != 1){
    resPos <- sapply(Seqs, function(x) words.pos(motif, x))
    resPos <- sapply(resPos, function(x) x+cenPos)
  }else{
    res <- sapply(Seqs, function(x) words.pos(motif, x))
    resPos <- list()
    if(is.null(dim(res))){
      resPos[[names(Seqs)]] <- res[1] + cenPos
    }else{
      resPos[[names(Seqs)]] <- res[,1] + cenPos
    }
  }
  resPos
}

#' @export
motifScan <- function(sequence, motif = "[ag][ag]ac[act]"){
  peakSeq <- read.fasta(file = sequence, as.string = T)
  motifPos <- searchMotifPos(sequence = sequence)
  kk <- unlist(lapply(X = motifPos, function(x) length(x)))
  motifPos <- motifPos[which(kk != 0)]
  
  ll <- 1
  resSeq <- list()
  for(i in 1:length(motifPos)){
    curID <- names(motifPos)[i]
    curPos <- motifPos[[i]]
    if(length(curPos) != 0){
      for(j in 1:length(curPos)){
        curSeq <- substr(peakSeq[[curID]], curPos[j]-2, curPos[j]+2)
        resSeq[[ll]] <- curSeq
        ll <- ll + 1
      }
    }
  }
  
  pfw <- matrix(0, 4, ncol = 5)
  rownames(pfw) <- c("a", "c", "g", "t")
  
  for(i in 1:ncol(pfw)){
    curSeq <- substr(resSeq, i, i)
    curFreq <- table(curSeq)/length(curSeq)
    curFreq <- curFreq[which(names(curFreq) != "N")]
    pfw[names(curFreq), i] <- curFreq
  }
  
  pwm <- makePWM(pfw)
  seqLogo(pwm)
  pwm
}

#' @export
motifDetect <- function(sequence, plot = T, ...){
  sequence <- readDNAStringSet(filepath = sequence, format = 'fasta')
  results <- GADEM(sequence, verbose = 1, ...)
  motifList <- results@motifList
  resList <- list()
  for(i in 1:length(motifList)){
    resList[[i]] <- motifList[[i]]@pwm
  }
  
  if(plot){
    for(i in 1:length(resList)){
      seqLogo(pwm = resList[[i]])
    }
    
  }
  
  resList
}


.getUTR <- function(curID, GTF){
  cur5UTRMat <- subset(GTF, GTF$V3 == "five_prime_utr" & GTF$V6 == curID)[,4:5]
  # cur5UTRMat <- as.numeric(GTF[which(GTF[,3] == 'five_prime_utr' & 
  #                                      GTF[,6] == curID), 4:5, drop = FALSE])
  cur3UTRMat <- subset(GTF, GTF$V3 == "three_prime_utr" & GTF$V6 == curID)[,4:5]
  curIDRange <- subset(GTF, GTF$V3 == "CDS" & GTF$V6 == curID)[,4:5]
  
  if(nrow(cur5UTRMat) >= 1){
    len.5 <- sum(cur5UTRMat$V5 - cur5UTRMat$V4 + 1)
  }else{
    len.5 <- 0
  }
  
  if(nrow(cur3UTRMat) >= 1){
    len.3 <- sum(cur3UTRMat$V5 - cur3UTRMat$V4 + 1)
  }else{
    len.3 <- 0
  }
  
  cds.len <- sum(curIDRange$V5 - curIDRange$V4 + 1)
  if(len.5 != 0){
    UTR5Range <- c(1, len.5)
    CDSRange <- c((len.5 + 1), (len.5 + cds.len))
  }else{
    UTR5Range <- c(0, 0)
    CDSRange <- c(1, cds.len)
  }
  
  if(len.3 != 0){
    UTR3Range <- c((cds.len + len.5 + 1), (len.5 + len.3 + cds.len))
  }else{
    UTR3Range <- c(0, 0)
  }
  res <- unlist(c(UTR5Range, CDSRange, UTR3Range))
}


#' @export
getUTR <- function(GTF, cpus = 1){
  GTF <- read.table(GTF, sep = "\t", header = F, 
                    stringsAsFactors = F, quote = "")
  GTF$V3[grep(pattern = "5|five", x = GTF$V3)] <- "five_prime_utr"
  GTF$V3[grep(pattern = "3|three", x = GTF$V3)] <- "three_prime_utr"
  GTF$V9 <- gsub(pattern = "\"", replacement = "", x = GTF$V9)
  
  curGTF <- subset(GTF, GTF$V3 == "CDS" | GTF$V3 == "five_prime_utr" | GTF$V3 == "three_prime_utr")
  curGTF$V6 <- sapply(curGTF$V9,  .extractTranscriptID)
  
  transcriptID <- unique(curGTF[,6])
  if(cpus == 1){
    UTRMat <- t(sapply(X = transcriptID, FUN = .getUTR, GTF = curGTF))
  }else{
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport(".getUTR", namespace = "PEA")
    UTRMat <- t(sfSapply(transcriptID, .getUTR, GTF = curGTF))
    sfStop()
  }
  
  UTRMat <- cbind(rownames(UTRMat), UTRMat)
  colnames(UTRMat) <- c("cDNA_ID",
                        "five_UTR_Start",
                        "five_UTR_End",
                        "CDS_Start",
                        "CDS_End",
                        "three_UTR_Start",
                        "three_UTR_End")
  rownames(UTRMat) <- UTRMat[,1]
  UTRMat <- UTRMat[, -1]
  UTRMat
}


.G2T <- function(posVec, exonGTF){
  curCHR <- posVec[1]
  curExon <- exonGTF[which(exonGTF[,1] == curCHR), ]
  curPos1 <- as.numeric(posVec[2])
  curPos2 <- as.numeric(posVec[3])
  
  index1 <- which((as.numeric(curExon[,4]) <= curPos1) & (as.numeric(curExon[,5]) >= curPos1))
  index2 <- which((as.numeric(curExon[,4]) <= curPos2) & (as.numeric(curExon[,5]) >= curPos2))
  index <- unique(c(index1, index2))
  
  if(length(index) == 0){
    resMat <- NULL
  }else{
    resSamples1 <- NULL
    resSamples2 <- NULL
    resMat <- matrix(NA, nrow = length(index), ncol = 3)
    for(j in 1:length(index)){
      curidx <- index[j]
      curExonParen <- curExon[curidx, 6]
      curStrand <- curExon[curidx, 7]
      curExonStart <- as.numeric(curExon[curidx, 4])
      curExonStop <- as.numeric(curExon[curidx, 5])
      curTranscript <- curExon[which(curExon[,6] == curExonParen), , drop = FALSE]
      if(curStrand == "-"){
        curTranscript <- curTranscript[order(as.numeric(curTranscript[,4]), decreasing = T), , drop = FALSE]
      }
      
      if(is.element(curidx, index1) & is.element(curidx, index2)){
        curExonidx1 <- which(as.numeric(curTranscript[,4]) <= curPos1 & as.numeric(curTranscript[,5]) >= curPos1)
        curExonidx2 <- which(as.numeric(curTranscript[,4]) <= curPos2 & as.numeric(curTranscript[,5]) >= curPos2)
        
        if(curStrand == "-"){
          tmpPos1 <- as.numeric(curTranscript[curExonidx1, 5]) - curPos2 + 1
          tmpPos2 <- as.numeric(curTranscript[curExonidx2, 5]) - curPos1 + 1
        }else{
          tmpPos1 <- curPos1 - as.numeric(curTranscript[curExonidx1, 4]) + 1
          tmpPos2 <- curPos2 - as.numeric(curTranscript[curExonidx2, 4]) + 1
        }
        
        ExonLen <- as.numeric(curTranscript[,5]) - as.numeric(curTranscript[,4]) + 1
        ExonLen <- c(0, ExonLen)
        resPos1 <- sum(ExonLen[1:curExonidx1]) + tmpPos1
        resPos2 <- sum(ExonLen[1:curExonidx2]) + tmpPos2
        # resSamples1 <- c(resSamples1, paste0(curExonParen, "_", as.numeric(resPos1)))
        # resSamples2 <- c(resSamples2, paste0(curExonParen, "_", as.numeric(resPos2)))
      }else if(is.element(curidx, index1)){
        curExonidx1 <- which(as.numeric(curTranscript[,4]) <= curPos1 & as.numeric(curTranscript[,5]) >= curPos1)
        curExonidx2 <- curExonidx1
        
        if(curStrand == "-"){
          tmpPos1 <- as.numeric(curTranscript[curExonidx1, 5]) - curPos1 + 1
          tmpPos2 <- 0
        }else{
          tmpPos1 <- curPos1 - as.numeric(curTranscript[curExonidx1, 4]) + 1
          tmpPos2 <- as.numeric(curTranscript[curExonidx2, 5]) - as.numeric(curTranscript[curExonidx2, 4]) + 1
        }
        
        ExonLen <- as.numeric(curTranscript[,5]) - as.numeric(curTranscript[,4]) + 1
        ExonLen <- c(0, ExonLen)
        resPos1 <- sum(ExonLen[1:curExonidx1]) + tmpPos1
        resPos2 <- sum(ExonLen[1:curExonidx2]) + tmpPos2
        # resSamples1 <- c(resSamples1, paste0(curExonParen, "_", as.numeric(resPos1)))
        # resSamples2 <- c(resSamples2, paste0(curExonParen, "_", as.numeric(resPos2)))
      }else{
        curExonidx2 <- which(as.numeric(curTranscript[,4]) <= curPos2 & as.numeric(curTranscript[,5]) >= curPos2)
        curExonidx1 <- curExonidx2
        
        if(curStrand == "-"){
          tmpPos1 <- as.numeric(curTranscript[curExonidx1, 5]) - as.numeric(curTranscript[curExonidx1, 4]) + 1
          tmpPos2 <- as.numeric(curTranscript[curExonidx2, 5]) - curPos2 + 1
        }else{
          tmpPos1 <- 0
          tmpPos2 <- curPos2 - as.numeric(curTranscript[curExonidx2, 4]) + 1
        }
        
        ExonLen <- as.numeric(curTranscript[,5]) - as.numeric(curTranscript[,4]) + 1
        ExonLen <- c(0, ExonLen)
        resPos1 <- sum(ExonLen[1:curExonidx1]) + tmpPos1
        resPos2 <- sum(ExonLen[1:curExonidx2]) + tmpPos2
        
      }
      
      resStart <- min(c(resPos1, resPos2))
      resStop <- max(c(resPos1, resPos2))
      resMat[j,] <- c(curExonParen, resStart, resStop)
      
    }
    
    
    
  }
  resMat
}

.extractTranscriptID <- function(x){
  res <- unlist(strsplit(x, ";"))
  idx <- grep("transcript_id", res)
  res <- res[idx]
  resID <- gsub(pattern = "transcript_id| ", replacement = "", x = res)
  resID
}

#' @export
G2T <- function(bedPos, GTF){
  GTF <- read.table(GTF, sep = "\t", quote = "", header = F, stringsAsFactors = F)
  interChr <- intersect(GTF[,1], bedPos[,1])
  if(length(interChr) == 0){
    cat("The chromosomes are not consistent in the GTF and bedPos, please provide right input!")
    return(NULL)
  }
  GTF$V9 <- gsub(pattern = "\"", replacement = "", x = GTF$V9)
  exonGTF <- GTF[which(GTF$V3 == "exon"),]
  exonGTF[,6] <- apply(exonGTF[, 9, drop = FALSE], 1, .extractTranscriptID)
  if(ncol(bedPos) > 3){
    bedPos <- bedPos[,1:3]
  }
  resPos <- apply(bedPos, 1, .G2T, exonGTF = exonGTF)
  resPos <- do.call(rbind, resPos)
  resPos
}



.geneID <- function(GTF){
  GTF <- read.table(file = GTF, header = F, sep = "\t", stringsAsFactors = F, quote = "")
  GTF$V9 <- gsub(pattern = "\"", replacement = "", x = GTF$V9)
  GTF <- subset(x = GTF, GTF$V3 == "transcript")
  curMat <- strsplit(GTF$V9, ";")
  curMat <- lapply(curMat, FUN = function(x) x[1:2])
  curMat <- do.call(what = rbind, curMat)
  curMat[,1] <- substr(curMat[,1], 9, nchar(curMat[,1])) 
  curMat[,2] <- substr(curMat[,2], 15, nchar(curMat[,2])) 
  resMat <- curMat[,1:2]
  resMat[,1] <- gsub(pattern = " ", replacement = "", resMat[,1])
  resMat[,2] <- gsub(pattern = " ", replacement = "", resMat[,2])
  resMat
}

.T2G <- function(transcript, geneMat){
  if(is.element(transcript, rownames(geneMat))){
    res <- geneMat[transcript, 1]
  }else{
    res <- NA
  }
  res
}

#' @export
CMRAnnotation <- function(cmrMat = NULL, genomic = F, UTRMat = NULL, GTF = NULL, SNR = T,
                          annotation = c("location", "motifScan", "motifDetect", "GO"),
                          cmrSeq = NULL, RNAseq = NULL, motifPos = NULL, plot = T, ...){
  
  if(length(annotation) > 1){
    cat("Warnings: multiple annotation was provided, the first one will be used!")
    annotation <- annotation[1]
  }
  
  
  if(genomic){
    cmrMat <- G2T(bedPos = cmrMat, GTF = GTF)
    #peakMat <- do.call(peakMat, rbind)
  }else{
    cmrMat <- cmrMat
  }
  
  if(annotation == "location"){
    if(is.null(GTF)){
      stop("Please provide the GTF!", "\n")
    }
    
    geneID <- .geneID(GTF = GTF)
    geneID <- unique(geneID)
    geneID[,2] <-  gsub(" ", "", x = geneID[,2], fixed = TRUE)
    rownames(geneID) <- geneID[,2]
    
    if(is.null(UTRMat)){
      UTRMat <- getUTR(GTF = GTF)
    }
    
    
    class(UTRMat) <- "numeric"
    
    if(!SNR){
      if(is.null(RNAseq)){
        stop("Please provide the RNA sequence!")
      }
      peaks <- cmrMat
      if(is.null(motifPos)){
        motifPos <- searchMotifPos(sequence = RNAseq, ...)
      }
      
      .convertSNR <- function(inputVec, motifPos){
        curID <- inputVec[1]
        curStart <- as.numeric(inputVec[2])
        curEnd <- as.numeric(inputVec[3])
        curMotif <- motifPos[[curID]]
        curPos <- curMotif[which(curMotif >= curStart & curMotif <= curEnd)]
        tmpMat <- matrix(NA, nrow = length(curPos), 2)
        tmpMat[,1] <- curID
        tmpMat[,2] <- curPos
        tmpMat
      }
      resPosMat <- apply(cmrMat, 1, .convertSNR, motifPos = motifPos)
      resPosMat <- do.call(rbind, resPosMat)
      colnames(resPosMat) <- c("IDs", "Position")
      peakMat <- resPosMat
    }else{
      peakMat <- cmrMat
    }
    
    ################CMR distribution################################
    .UTRPosition <- function(inputVec, UTRMat){
      posSampleID <- inputVec[1]
      posSamplePos <- as.numeric(inputVec[2])
      curUTR <- as.numeric(UTRMat[posSampleID, ])
      if(posSamplePos >= curUTR[5] & posSamplePos <= curUTR[6]){
        res <- "three prime UTR"
      }else if(posSamplePos >= curUTR[3] & posSamplePos <= curUTR[4]){
        res <- "CDS"
      }else{
        res <- "five prime UTR"
      }
      res
    }
    peakMat <- peakMat[which(peakMat[,1] %in% rownames(UTRMat)),]
    finalPosition <- apply(peakMat, 1, .UTRPosition, UTRMat = UTRMat)
    
    peakMat <- cbind(peakMat, finalPosition)
    resPos <- table(finalPosition)
    
    ########################CMR normalized distribution################
    .normalizeDis <- function(inputVec, UTRMat){
      curID <- inputVec[1]
      curPos <- as.numeric(inputVec[2])
      curRegion <- inputVec[3]
      
      if(curRegion == "five prime UTR"){
        curRes <- curPos/UTRMat[curID, 2]
      }else if(curRegion == "CDS"){
        curRes <- (curPos - UTRMat[curID,2] + 1)/((UTRMat[curID, 4] - UTRMat[curID, 3]) + 1) + 1
      }else{
        curRes <- (curPos - UTRMat[curID,5] + 1)/((UTRMat[curID, 6] - UTRMat[curID, 5]) + 1) + 2
      }
      curRes
    }
    positive.sample.position <- apply(peakMat, 1, .normalizeDis, UTRMat = UTRMat)
    
    if(plot){
      if(SNR){
        par(mfrow = c(1,3))
        par(mar=c(2, 2, 2, 2))   
        pie(resPos, col = c('yellow', "green", 'red'), labels = names(resPos),
            main = "CMR distribution along the transcript")
        
        plot(density(positive.sample.position), main = 
               "Distribution of CMR in the cDNA",
             col = "red", lwd = 2, xaxt = "n")
        rug(seq(0, 1, 0.001), col = "lightgreen")
        rug(seq(1, 2, 0.001), col = "cadetblue3")
        rug(seq(2, 3, 0.001), col = "red")
        legend("topleft", col = c("lightgreen", "cadetblue3", "red"), lwd = c(5, 5, 5), 
               legend = c("5'UTR", "CDS", "3'UTR"))
        aa <- table(peakMat[,1])
        bb <- table(aa)
        barplot(bb, col = rainbow(length(bb)),
                main = "Transcripts with different CMRs.")
      }else{
        par(mfrow = c(2,2))
        par(mar=c(2, 2, 2, 2))   
        pie(resPos, col = c('yellow', "green", 'red'), labels = names(resPos),
            main = "CMR distribution along the transcript")
        
        plot(density(positive.sample.position), main = 
               "Distribution of CMR in the cDNA",
             col = "red", lwd = 2, xaxt = "n")
        rug(seq(0, 1, 0.001), col = "lightgreen")
        rug(seq(1, 2, 0.001), col = "cadetblue3")
        rug(seq(2, 3, 0.001), col = "red")
        legend("topleft", col = c("lightgreen", "cadetblue3", "red"), lwd = c(5, 5, 5), 
               legend = c("5'UTR", "CDS", "3'UTR"))
        .peakmotifNumber <- function(inputVec, motifPos){
          curID <- inputVec[1]
          curMotif <- motifPos[[curID]]
          motifNumber <- length(which(curMotif >= as.numeric(inputVec[2]) & curMotif <= as.numeric(inputVec[3])))
          motifNumber
        }
        peakMotifNumber <- apply(cmrMat, 1, .peakmotifNumber, motifPos = motifPos)
        tt <- table(peakMotifNumber)
        barplot(tt, col = rainbow(length(tt)), main = "Peaks with different CMRs")
        aa <- table(peakMat[,1])
        bb <- table(aa)
        barplot(bb, col = rainbow(length(bb)), main = "Transcripts with different CMRs.")
      }
    }
    
    resList <- list(cmrMat = cmrMat, position = finalPosition,
                    distribution = positive.sample.position)
    return(resList)
  }
  
  
  if(annotation == "motifScan"){
    if(is.null(cmrSeq)){
      stop("Please provide the CMR-related sequences.")
    }
    
    results <- motifScan(sequence = cmrSeq, ...)
    return(results)
  }
  
  if(annotation == "motifDetect"){
    if(is.null(cmrSeq)){
      stop("Please provide the CMR-related sequences.")
    }
    
    results <- motifDetect(sequence = cmrSeq, ...)
    return(results)
  }
  
  if(annotation == "GO"){
    if(is.null(GTF)){
      stop("Please provide the GTF file name!")
    }
    geneID <- .geneID(GTF = GTF)
    resGene <- geneID[which(geneID[,2] %in% cmrMat[,1]),1]
    resGO <- runTopGO(geneID = resGene, ...)
    return(resGO)
  }
  resList
}

