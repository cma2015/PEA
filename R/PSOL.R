############################################################
##general functions in ML


#Function: get sample index for cv cross validation
.cvSampleIndex <- function( len, cross = 5, seed = 1 ) {
  
  cv <- cross
  sample_matrix <- matrix(0, nrow = len, ncol = cv)
  colnames(sample_matrix) <- paste("cv", c(1:cv), sep = "" )
  
  #random samples 
  set.seed(seed)
  index <- sample(1:len, len, replace = FALSE )
  step = floor( len/cv )
  
  start <- NULL
  end <- NULL
  train_lens <- rep(0, cv)
  for( i in c(1:cv) ) {
    start <- step*(i-1) + 1
    end <- start + step - 1
    if( i == cv ) 
      end <- len
    
    train <- index[-c(start:end)]
    test <- index[start:end]
    train_lens[i] <- length(train)
    
    sample_matrix[,i] <- c(train, test)
  }#end for i
  
  return( list( train_lens = train_lens, sample_matrix = sample_matrix))
}



#' @export
classifier <- function( method = c("randomForest", "svm"), featureMat, positiveSamples, 
                        negativeSamples, ...) {
  call <- match.call()
  
  if( length(method) > 1){
    method <- method[1]
  } 
    
  
  if( is.null(rownames(featureMat) ) )
    stop("Error: no row names (i.e., sample IDs) were assigned for featureMat." )
  if( is.null(colnames(featureMat) ) )
    stop("Error: no colnames were defined for featureMat." )
  
  positiveSamples <- intersect( rownames(featureMat), positiveSamples )
  negativeSamples <- intersect( rownames(featureMat), negativeSamples )
  posLen <- length(positiveSamples)
  negLen <- length(negativeSamples)
  if( posLen == 0 )
    stop("Error: no positive samples included in featureMat." )
  if( negLen == 0 )
    stop("Error: no negative samples were included in featureMat." )
  
  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(positiveSamples, negativeSamples), ] )
  tmpData <- cbind( fmat, label )
  colnames(tmpData) <- c(colnames(fmat), "Class")
  if( method == "randomForest" ) {
    obj <- randomForest(x = fmat, y = factor(label), ... )
  }else{
    obj <- svm(x = fmat, y = factor(label), ... )
  }
  obj
}






.predictor <- function( method = c("randomForest", "svm"), classifier, featureMat ) {
  
  if(length(method) > 1){
    method <- method[1]
  }
  
  if( method == "randomForest") {
    res <- predict(classifier, data.frame(featureMat), type= "vote" )[,"1"]
  }else {
    res <- predict( classifier, data.frame(featureMat), type = "raw") 
  }
  names(res) <- rownames(featureMat)
  res
}




.find_ClassifierWithMaxAUC <- function( cvRes ) {
  
  classifier <- NA
  maxAUC <- 0
  for( i in 1:length(cvRes) ) {
    res <- cvRes[[i]]
    if( res$test.AUC > maxAUC) {
      maxAUC <- res$test.AUC
      classifier <- res$classifier
    }
  }#end for i
  
  return( list(maxAUC = maxAUC, classifier = classifier))
}



.obtain_CV_AUCMat <- function( cvRes ) {
  cv <- length(cvRes)
  AUCMat <- matrix(0, nrow = cv, ncol = 2 )
  rownames(AUCMat) <- paste( "cv", 1:cv, sep = "" )
  colnames(AUCMat) <- c("trainingSet", "testingSet")
  
  for( i in 1:cv ) {
    res <- cvRes[[i]]
    AUCMat[i,2] <- res$test.AUC
    AUCMat[i,1] <- res$train.AUC
  }#end for i
  
  AUCMat
}

##get system time for seed and then generate random index
.randomSeed <- function() {
  curtime <- format(Sys.time(), "%H:%M:%OS4")
  XXX <- unlist(strsplit(curtime, ":"))
  curtimeidx <- (as.numeric(XXX[1])*3600 + as.numeric(XXX[2])*60 + as.numeric(XXX[3]))*10000
  curtimeidx
}



.one_cross_validation <- function( cv, method, featureMat, positives, negatives, posSample_cv, negSample_cv, balanced = TRUE, ratio = 10, ... ) {
  call <- match.call()
  j <- cv
  
  #for train samples
  train_genes_p <- positives[ (posSample_cv$sample_matrix[,j][1:posSample_cv$train_lens[j]] ) ]
  test_genes_p <- positives[ (posSample_cv$sample_matrix[,j][-c(1:posSample_cv$train_lens[j])]) ]
  
  #trained negatives randomly selected, and tested on all negatives
  train_genes_n <- negatives[(negSample_cv$sample_matrix[,j][1:negSample_cv$train_lens[j]] ) ]
  test_genes_n <- negatives[ (negSample_cv$sample_matrix[,j][-c(1:negSample_cv$train_lens[j])]) ]
  
  #select part of train_genes_n
  if( balanced == TRUE ) {
    if( length(train_genes_n) > ratio*length(train_genes_p) ) {
      train_genes_n <- train_genes_n[sample(1:length(train_genes_n), replace=FALSE)[1:(ratio*length(train_genes_p))]]
    }
  }
  
  
  
  obj <- classifier( method = method, featureMat = featureMat, positiveSamples = train_genes_p, negativeSamples = train_genes_n, ... )
  bestmodel <- obj
  
  positives.train.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[train_genes_p,])
  negatives.train.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[train_genes_n,])
  positives.test.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[test_genes_p,])
  negatives.test.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[test_genes_n,])
  
  
  
  train.AUC <- roc( c(rep(1, length(train_genes_p)), rep(0, length(train_genes_n))), 
                    c(positives.train.score, negatives.train.score) )$auc[1]
  test.AUC <- roc( c(rep(1, length(test_genes_p)), rep(0, length(test_genes_n))), 
                   c(positives.test.score, negatives.test.score) )$auc[1]
  
  res <- ( list( positves.train = train_genes_p, negatives.train = train_genes_n, 
                 positives.test = test_genes_p, negatives.test = test_genes_n, 
                 ml = method, classifier = bestmodel, 
                 positives.train.score = positives.train.score,
                 negatives.train.score = negatives.train.score,
                 positives.test.score = positives.test.score,
                 negatives.test.score = negatives.test.score,
                 train.AUC = train.AUC,
                 test.AUC = test.AUC) )
  
  res
}



#' @export
cross_validation <- function( seed = 1, method = c("randomForest", "svm"), 
                              featureMat, positives, negatives, cross = 5, 
                              cpus = 1, ... ){
  
  call <- match.call()
  
  if( length(method) > 1){
    method <- method[1]
  } 
  
  #sample index for cv
  posSample_cv <- .cvSampleIndex(length(positives), cross = cross, seed = seed)
  negSample_cv <- .cvSampleIndex(length(negatives), cross = cross, seed = seed)
  
  cvRes <- list()
  if( cpus > 1 ) {
    #require(snowfall)
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport("classifier", namespace = "PEA")
    sfExport(".predictor", namespace = "PEA")
    sfExport(".one_cross_validation", namespace = "PEA")
    sfLibrary( "pROC", character.only = TRUE)
    sfLibrary( "e1071", character.only = TRUE)
    sfLibrary( "randomForest", character.only = TRUE )
    
    cvRes <- sfApply( matrix(1:cross, ncol = 1), 1,  .one_cross_validation, method = method, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv, ...)
    sfStop()
  }else {
    for( j in 1:cross ) {
      cvRes[[j]] <- .one_cross_validation( cv = j, method = method, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv, ... )
    }
  }
  cvRes
}


.disMat <- function(inputVec, featureMatrix){
  curRes <- dist(featureMatrix[inputVec, ])
  curRes
}


.calDisMat <- function(speSample, featureMatrix, resDic){
  tmpMat <- matrix(NA, nrow = nrow(featureMatrix), ncol = 2)
  tmpMat[,1] <- 1:nrow(featureMatrix)
  tmpMat[,2] <- speSample
  resVec <- apply(tmpMat, 1, .disMat, featureMatrix)
  fileName <- paste0(rownames(featureMatrix)[speSample], ".RData")
  save(resVec, file = paste0(resDic, fileName))
}


.calculateDistance <- function(sampLength, featureMatrix, cpus, resDic){
  if(cpus > 1){
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport(".calDisMat", namespace = "PEA")
    sfExport(".disMat", namespace = "PEA")
    results <- sfApply(matrix(1:sampLength, ncol = 1), 1, .calDisMat, featureMatrix =  featureMatrix, resDic = resDic)
  }else{
    for(i in 1:sampLength){
      results <- .calDisMat(speSample = i, featureMatrix = featureMatrix, resDic = resDic)
    }
  }
}


#According to the paper Bioinformatics, 2006, 22:21 2590-2596, maximum distance minimum redundancy negative set is selected
#' @export
.PSOL_InitialNegativeSelection <- function( featureMatrix,  positives, unlabels, 
                                           negNum = length(positives), cpus = 1,
                                           PSOLResDic) {

   
  if( negNum > length(unlabels) )
    stop("Error: negNum is larger than the number of unlables")
  if( is.null( rownames(featureMatrix)) )
    stop("Error: no rownames for featureMatrix")

  ##create PSOLResDic
  
  resDic <- paste0(PSOLResDic, "adjMat/")
  dir.create( path = resDic, showWarnings = FALSE)
  results <- .calculateDistance(sampLength = nrow(featureMatrix), 
                                featureMatrix = featureMatrix, cpus = cpus, 
                                resDic = resDic)
  
  
  sampleNum <- nrow(featureMatrix)
  sampleName <- rownames(featureMatrix)
  
  ###############Creating bigmatrix##########################################
  cat("Start creating bigmatrix......")
  adjmat <- big.matrix(sampleNum, sampleNum, type="double", init=0.0, dimnames=list(sampleName, sampleName),
                         backingfile = "adjmat_backingfile", backingpath = PSOLResDic, descriptorfile = "adjmat_descriptfile")
  
  pb <- txtProgressBar(min = 0, max = length(sampleName), style = 3,width = 75)
  for(i in 1:length(sampleName)){
    curDir <- paste0(resDic, sampleName[i], ".RData")
    load(curDir)
    adjmat[i,] <- resVec
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  negative_new <- NULL  
  
  unlabel_new <- unlabels
  min_up <- apply( adjmat[unlabel_new, positives], 1, min )
  names(min_up) <- unlabel_new
  
  cat( "negativeSampleSize:")
  negCount <- 0
  while( negCount <= negNum ) {
    
    if( negCount%%10 == 0 )
      cat(negCount, "..")
    
    index <- -1

    #start, no negative selected
    if( negCount == 0 ) { 
            
      #find one unlabel sample with max EDistance to positives
      index <- which( min_up == max(min_up) )[1]      
    }else {
      
      #comput E distance between unlabels and positives, unlables and negatives
      res_un <- adjmat[unlabel_new, negative_new ]
      if( negCount == 1 )
        res_un <- matrix( res_un, ncol = 1 )
     
      #maximum distance minimum redudancy
      sum_un <- apply(res_un, 1, sum)
      res <- t(min_up)*sum_un
      index <- which( res == max(res) )[1]
    }

    #record selected unlabel samples as negatives
    negCount <- negCount + 1
    negative_new <- c( negative_new, unlabel_new[index])
      
    #remove recorded negatives from featureMatrix.unlabel
    unlabel_new <- unlabel_new[-index]
    min_up <- min_up[-index]
    
    if( negCount == negNum )
      break
    
  }
  cat("\n")

  res <- list( positives = positives, negatives = negative_new, unlabels = unlabel_new)
  save( res, file = paste( PSOLResDic, "PSOL_InitialNegativeSelection_Res.RData", sep = "" ) )
  res
}





####PSOL Negative Exapansion
#' @export
.PSOL_NegativeExpansion <- function( featureMat, positives, negatives, unlabels,
                                     cpus = 1,  iterator = 50, cross = 5, 
                                     TPR = 0.98, method = c("randomForest", "svm"), 
                                     plot = TRUE, trace = TRUE, PSOLResDic, ... ) {

  call <- match.call()
  if( length(method) > 1){
    method <- method[1]
  } 
  
  ##create PSOLResDic
  dir.create( path = PSOLResDic, showWarnings = FALSE)

  if( TPR > 1.0 | TPR < 0 )
     stop("Error: TPR is a value ranged from 0 to 1.")
  if( is.null(rownames(featureMat)) )
    stop("Error: rownames should be given to featureMat")

  ##get the thresholdIdx-th top prediction score of positive samples as the threshold cutoff
  thresholdIdx <- floor( (1.0 - TPR)*length(positives) )
  if( thresholdIdx <= 0 )
    thresholdIdx <- 1

  finalUnlabels <- unlabels
  finalNegatives <- negatives
  negCount <- length(finalNegatives)
  
  numMat <- matrix(0, nrow = iterator, ncol = 5 )
  rownames(numMat) <- paste("Iter", 1:iterator, sep = "" )
  colnames(numMat) <- c("IterNo", "AUC_On_TrainingDataSet", "AUC_On_TestingDataSet", "Negative_Sample_Num", "Unlabeled_Sample_Num")
  numMat[,1] <- 1:iterator

  ##check samples
  if( length( setdiff( c( positives, negatives, unlabels), rownames(featureMat) ) ) > 0 ) {
    stop("Error: some samples not included in the featureMat.\n")
  }
  featureMat <- featureMat[c( positives, negatives, unlabels), ]
  
  
  zeroNumCount = 0
  iter <- 0
  while( negCount <= length(unlabels) ){
    
    iter <- iter + 1
    if( iter > iterator ){
      iter <- iter - 1
      break
    }
    #cross validation   
    permutRes <- cross_validation( seed = .randomSeed(), method = method, featureMat = featureMat, positives = positives, negatives = finalNegatives, cross = cross, cpus = cpus, ... )
    
    ##find classifier with max AUC, and re-calculate prediction scores of positive and unlable samples
    maxAUC_Classifer <- .find_ClassifierWithMaxAUC( permutRes )
    prediction.score <- .predictor( method = method, classifier = maxAUC_Classifer$classifier, featureMat = featureMat ) 
    positives.score <- sort( prediction.score[positives], decreasing = FALSE )
    negatives.score <- prediction.score[finalNegatives] 
    unlabels.score <- prediction.score[finalUnlabels]

    positive.Threshold <- as.numeric( positives.score[thresholdIdx] )  #threshold
    num <- length( which( unlabels.score < positive.Threshold) ) 
    finalUnlabels <- names(unlabels.score[which(unlabels.score > positive.Threshold)])
    finalNegatives <- unique( c( names(unlabels.score[which(unlabels.score <= positive.Threshold)]), finalNegatives) )  

    AUCMat <- .obtain_CV_AUCMat( permutRes )
    numMat[iter, 2] <- mean(AUCMat[,1]) #AUC on training dataset
    numMat[iter, 3] <- mean(AUCMat[,2]) #AUC on testing dataset
    numMat[iter, 4] <- length(negatives.score)
    numMat[iter, 5] <- length(unlabels.score)
    
    #plot density
#    if( plot == TRUE ) {
#       pdf( paste( PSOLResDic, "PSOL_Iteration_", iter, ".pdf", sep = ""), height= 5, width= 10)
#       par(mfrow=c(1,2))
#       density.p <- density(positives.score)
#       density.n <- density(unlabels.score)
#       xrange = range(density.p$x, density.n$x)
#       yrange = range(density.p$y, density.n$y)
#       plot( density.p, xlim = xrange, ylim = yrange, col = "red", xlab="prediction score", ylab = "density", main = paste("iterator times:", iter, sep="" ) )
#       lines(density.n, col = "black")
#       abline( v = positive.Threshold, col = "gray" )
    
       #plot differences of AUC
#       boxplot(AUCMat, ylim = c(0,1) )
#       dev.off()
#    } 
        
    #save result
    if( trace == TRUE ) {
       resultDir <- paste( PSOLResDic, "PSOL_Iteration_", iter, ".RData", sep = "")
       iterRes <- list( permutRes = permutRes, method = method, classifier = maxAUC_Classifer, 
                     predictionScores = prediction.score, negativeScores = negatives.score, 
                     unlabelScores = unlabels.score, threshold = positive.Threshold, 
                     positives = positives, negatives = names(negatives.score), 
                       unlabels = names(unlabels.score), #here negatives mean the "negatives" started at this iteration time. 
                     finalNegatives = finalNegatives, finalUnlabels = finalUnlabels )
       save( iterRes, file = resultDir )
    }
    
    cat( "\nPSOL_Iteration: ", iter, "\tAUC: ", numMat[iter, 3], "\tCurrentPosNum:", length(positives),  "\tCurrentNegNum: ", numMat[iter, 4], "\tCurrentUnlabelNum: ", numMat[iter, 5], "\tIncreased negatives Num: ", num, "\n")
     
  }##end while
  
  ##plot number distribution
  numMat <- numMat[1:iter,]
  write.table( numMat, paste( PSOLResDic, "PSOL_NegativeIncreasement.txt", sep="" ), sep = "\t", quote = F )
  if( plot == TRUE ) {
    

    pdf( paste( PSOLResDic, "PSOL_NegativeIncreasement.pdf", sep="" ), height= 10, width = 10 )  
    par(mar=c(5, 12, 4, 4) + 0.1)
    plot(numMat[,1], numMat[,3], axes=F, ylim=c(0,1.0), xlab="", ylab="",type="l",col="red", main="")
    points(numMat[,1],numMat[,3],pch = 20,col = "red", cex = 0.8)
    axis(2, ylim=c(0,1.0),col="red",lwd=2)
    mtext(2,text="AUC",line=2)
  
    par(new=T)
    plot(numMat[,1], numMat[,4], axes=F, ylim=c(0,max(numMat[,4])), xlab="", ylab="", type="l", col = "black", lty=2, main="",lwd=2)
    axis( 2, ylim=c(0,max(numMat[,4])), lwd=2, line=3.5, col = "black" )
    points(numMat[,1], numMat[,4], pch=20, col = "black", cex = 0.8)
    mtext(2,text="Number of \"filtered-out\" genes ",line=5.5)
  
    axis(1,numMat[,1] )
    mtext("Iteration Number",side=1,col="black",line=2)
    #legend("bottomright",legend=c("negatives","benchmark genes"),lty=c(1,2), col = c("black", "red") )
    dev.off()
  }
   
}

