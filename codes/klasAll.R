klasAll <- function(data1, L)
{
  require(MLInterfaces)
  mis_classified <- list()
  mis_classified[['knn']] <- NULL
  mis_classified[['svmI']] <- NULL
  mis_classified[['nnet']] <- NULL
  mis_classified[['randomForest']] <- NULL
  
  dane<-data.frame(data1,L)
  #knn
  for(i in 1:10)
  {
    klasyfikator = MLearn(L~., data=dane, knnI(k=3,l=0),xvalSpec("LOG", 5, balKfold.xvspec(5)))
    a<-confuMat(klasyfikator)
    good_predicted<-sum(diag(a))
    mis_classified[['knn']]<-c(mis_classified[['knn']], dim(dane)[1]-good_predicted)
    
    #svm
    klasyfikator = MLearn(L~., data=dane, svmI, xvalSpec("LOG", 5, balKfold.xvspec(5)))
    a<-confuMat(klasyfikator)
    good_predicted<-sum(diag(a))
    mis_classified[['svmI']]<-c(mis_classified[['svmI']],dim(dane)[1]-good_predicted)
    
    #nnet
    klasyfikator = MLearn(L~., data=dane, nnetI, xvalSpec("LOG", 5, balKfold.xvspec(5)),size=3, decay=.5, trace=FALSE )
    a<-confuMat(klasyfikator)
    good_predicted<-sum(diag(a))
    mis_classified[['nnet']]<- c(mis_classified[['nnet']],dim(dane)[1]-good_predicted)
    
    #random forest
    klasyfikator = MLearn(L~., data=dane, randomForestI, xvalSpec("LOG", 5, balKfold.xvspec(5)))
    a<-confuMat(klasyfikator)
    good_predicted<-sum(diag(a))
    mis_classified[['randomForest']]<- c(mis_classified[['randomForest']], dim(dane)[1]-good_predicted)
  }
  mis_classified
}

klasAll1 <- function(x, y)
{
  require(caret)
  k <- ifelse(length(y) < 20, 5, 10)
  if( length(y) < 10) k <- 3
  foldy <- createFolds(y, k = k)
  nr.sam <- length(unique(y))
  out <- data.frame(true.y = as.numeric(as.character(y)), yhat = NA)
  prob <- matrix(NA,length(y), 3)
  colnames(prob) <- c("X0", "X1", "X2")
  M <- 100
  for (g in 1:length(foldy))
  {
    x1 <- x[-foldy[[g]], ]
    y1 <- y[-foldy[[g]]]
    ens <- ensembleClassifier(x1, y1, M = M, algorithms=c("nnet", "rf", "rpart", "svm"), verbose = F)
  
    trueClass <- y[foldy[[g]]]
  
    pred <- predictEns(ens, x[foldy[[g]], ], trueClass)
    out$yhat[foldy[[g]]] <- as.numeric(as.character(pred$yhat))
    for (i in 1:nr.sam) prob[foldy[[g]],i] <- apply(pred$pred, 1, function(x) length(which(x == levels(y)[i]))/M)
  }
  out <- cbind(out, prob)
  out
}
