###########CMR (Chemical Modification of RNA) Calling
##########author: Jingjing Zhai, Chuang Ma
##########Contact: zhaijingjing603@gmail.com

.generateIndex <- function(alignment = c("bowtie", "bowtie2", "tophat",
                                         "tophat2", "hisat", "hisat2"), 
                           refGenome = NULL, softwareDir = NULL){
  
  if(is.null(refGenome)){
    stop("Please provide the reference genome!")
  }
  
  if(alignment == "tophat" | alignment == 'bowtie'){
    tmpPara <- paste0("bowtie-build ", refGenome, " ", refGenome)
  }else if(alignment == "tophat2" | alignment == 'bowtie2'){
    tmpPara <- paste0("bowtie2-build ", refGenome, " ", refGenome)
  }else{
    tmpPara <- paste0(alignment, "-build ", refGenome, " ", refGenome)
  }
  
  tmpPara
}

#' @export
readMapping <- function(alignment = c("bowtie", "bowtie2", "tophat", 
                                    "tophat2", "hisat", "hisat2"), 
                          fq, index = NULL, refGenome = NULL,
                          paired = FALSE, softwareDir = NULL, ...){
  
  if(length(alignment) > 1){
    alignment <- alignment[1]
    cat("Multiple alignmen methods are provided, only the first one will be used!\n")
  }
  
  
  if(!is.null(softwareDir)){
    alignment <- paste0(softwareDir, alignment)
  }
  
  inputPara <- list(...)
  if(length(inputPara) != 0){
    alignPara <- inputPara[[1]]
  }else{
    alignPara <- NULL
  }
  
  
  cat("Start using ", alignment, " to align reads to reference genome......\n")

  if(is.null(index)){
    tmpPara <- .generateIndex(alignment = alignment, refGenome = refGenome, 
                              softwareDir = softwareDir)
    system(tmpPara)
    index <- refGenome
  }
  

  if(!paired){
    #SRAName <- unlist(strsplit(x = fq[1], "/"))
    #SRAName <- SRAName[length(SRAName)]
    #SRAName <- unlist(strsplit(x = SRAName, split = ".", fixed = TRUE))[1]
    SRAName <- unlist(strsplit(fq, ".", fixed = T))[1]
    ######bowtie alignment######
    if(alignment == 'bowtie'){
      alignmentPara <- paste0(alignment, " ", alignPara, " ", index, " ", fq,
                              " ", SRAName, "_bowtie.sam")
      resDic <- paste0(SRAName, "_bowtie.sam")
    }
    
    ######bowtie2 alignment######
    if(alignment == 'bowtie2'){
      alignmentPara <- paste0(alignment, " ", alignPara, " -x ", index, 
                              " -U ", fq, " -S ", SRAName, "_bowtie.sam")
      resDic <- paste0(SRAName, "_bowtie2.sam")
    }
    
    ######tophat alignment######
    if(alignment == 'tophat'){
      alignmentPara <- paste0(alignment, " ", alignPara, " --bowtie1 ", " -o ",
                              SRAName, "_tophat ", index, " ", fq)
      resDic <- paste0(SRAName, "_tophat/accepted_hits.bam")
    }
    
    ######tophat2 alignment######
    if(alignment == 'tophat2'){
      alignmentPara <- paste0(alignment, " ", alignPara, " -o ", SRAName, 
                              "_tophat2 ", index, " ", fq[1], " ", fq[2])
      resDic <- paste0(SRAName, "_tophat2/accepted_hits.bam")
    }
    
    if(alignment == 'hisat' | alignment == "hisat2"){
      alignmentPara <- paste0(alignment, " ", alignPara, " -x ", index, 
                              " -U ", fq, " -S ", SRAName, "_", alignment, ".sam")
      resDic <- paste0(SRAName, "_", alignment, ".sam")
    }
    
  }else{
    ####Two fastq files should be provided for paired alignment
    if(length(fq) != 2){
      stop("Paired-end alignment should provide two fastq files!")
    }
    
    SRAName <- unlist(strsplit(fq[1], split = "_", fixed = T))[1]
    ######bowtie alignment######
    if(alignment == 'bowtie'){
      alignmentPara <- paste0(alignment, " ", alignPara, " ", index, " -1 ",
                              fq[1], " -2 ", fq[2], " ", SRAName, "_bowtie.sam")
      resDic <- paste0(SRAName, "_bowtie.sam")
    }
    
    ######bowtie2 alignment######
    if(alignment == 'bowtie2'){
      alignmentPara <- paste0(alignment, " ", alignPara, " -x ", index, " -1 ", 
                              fq[1], " -2 ", fq[2], " -S ", SRAName, "_bowtie2.sam")
      resDic <- paste0(SRAName, "_bowtie2.sam")
    }
    
    ######tophat alignment######
    if(alignment == 'tophat'){
      alignmentPara <- paste0(alignment, " ", alignPara, "  --bowtie1", " -o ",
                              SRAName, "_tophat ", index, " ", fq[1], " ", fq[2])
      resDic <- paste0(SRAName, "_tophat/accepted_hits.bam")
    }
    
    ######tophat2 alignment######
    if(alignment == 'tophat2'){
      alignmentPara <- paste0(alignment, " ", alignPara, " -o ", SRAName, 
                              "_tophat2 ", index, " ", fq[1], " ", fq[2])
      resDic <- paste0(SRAName, "_tophat2/accepted_hits.bam")
    }
    
    if(alignment == 'hisat' | alignment == "hisat2"){
      alignmentPara <- paste0(alignment, " ", alignPara, " -x ", index, " -1 ", fq[1],
                              " -2 ", fq[2], " -S ", SRAName, "_", alignment, ".sam")
      resDic <- paste0(SRAName, "_", alignment, ".sam")
    }
  }
  system(alignmentPara)
  resList <- list()
  resList[['command']] <- alignmentPara
  resList[["resDir"]] <- resDic
  resList
}


#' @export
QCBAM <- function(BAM, quality = c(35, 30), paired = F){
  
  if(paired){
    mapqFilter <- quality[2]
  }else{
    mapqFilter <- quality[1]
  }
  
  p1 <- ScanBamParam(mapqFilter = mapqFilter, what = scanBamWhat())
  res <- scanBam(file = BAM, param = p1)
  res
}
