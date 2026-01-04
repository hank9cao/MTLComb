rm(list=ls())
gc()


#######################################################################################################################
#
#                             Cross-cohort prediction: Train on cohort 1 and test on cohort 2
#
#######################################################################################################################



dat_cohort1 = readRDS("./0_data_cohort1.rds")
dat_cohort2 = readRDS("./0_data_cohort2.rds")

out_cohort1 = readRDS("./0_out_cohort1.rds")
out_cohort2 = readRDS("./0_out_cohort2.rds")


source('../../MTLComb/MTL_L21.R')
source('../../MTLComb/MTLComb.R')
library(glmnet)
library(randomForest)
library(e1071)
library(doParallel)
library(foreach)
registerDoParallel(10)


rData = foreach (rep = 1:10)%dopar%{
  print(rep)
  
  predAucVec = vector()
  coefMat = vector()
  
  ################################################################################
  # Machine learning analysis on all variables
  ################################################################################
  X1=dat_cohort1
  X1$c_dialyse=NULL
  X1$event=NULL
  X1=as.matrix(X1)
  X1=apply(X1, 2, function(x)scale(x))
  Y1=dat_cohort1$event
  
  X2=dat_cohort2
  X2$c_dialyse=NULL
  X2$event=NULL
  X2=as.matrix(X2)
  X2 = apply(X2, 2, function(x)scale(x))
  Y2=dat_cohort2$event
  
  cvR=cv.glmnet(x=X1, y=Y1, nfolds = 10, family="binomial", alpha = 0)
  fitRidge = glmnet(x=X1, y=Y1, family="binomial", lambda = cvR$lambda.min, alpha = 0)
  predScore=predict(fitRidge, X2)
  predAucVec = c(predAucVec, as.numeric(pROC::auc(as.factor(Y2), predScore)))
  
  cvR=cv.glmnet(x=X1, y=Y1, nfolds = 10, family="binomial", alpha = 1)
  fitLasso = glmnet(x=X1, y=Y1, family="binomial", lambda = cvR$lambda.min, alpha = 1)
  predScoreLasso=predict(fitLasso, X2)
  predAucVec = c(predAucVec, as.numeric(pROC::auc(as.factor(Y2), predScoreLasso)))
  
  fitRF = randomForest(y = as.factor(Y1), x = X1, ntree=5000)
  predScoreRF = predict(fitRF, X2, type="prob")[,2]
  predAucVec = c(predAucVec, as.numeric(pROC::auc(as.factor(Y2), predScoreRF)))
  
  fitSVM = svm(y=as.factor(Y1), x=X1, kernel = "linear", cost = 1, scale = FALSE, probability=TRUE)
  predScoreSVM = attr(predict(fitSVM, X2,probability=TRUE), "probabilities")[,2]
  predAucVec = c(predAucVec, as.numeric(pROC::auc(as.factor(Y2), predScoreSVM)))
  
  fitSVM2 = svm(y=as.factor(Y1), x=X1, kernel = "radial", cost = 1, scale = FALSE, probability=TRUE)
  predScoreSVM2 = attr(predict(fitSVM2, X2,probability=TRUE), "probabilities")[,2]
  predAucVec = c(predAucVec, as.numeric(pROC::auc(as.factor(Y2), predScoreSVM2)))
  
  fitSVM3 = svm(y=as.factor(Y1), x=X1, kernel = "polynomial", cost = 1, scale = FALSE, probability=TRUE)
  predScoreSVM3 = attr(predict(fitSVM3, X2,probability=TRUE), "probabilities")[,2]
  predAucVec = c(predAucVec, as.numeric(pROC::auc(as.factor(Y2), predScoreSVM3)))
  
  
  ################################################################################
  # MTLComb-feat analysis with inductive CV
  ################################################################################
  X1=dat_cohort1
  X1$c_dialyse=NULL
  X1$event=NULL
  X1=as.matrix(X1)
  xx=apply(X1, 2, function(x)scale(x))
  X=list(X1=xx, X2=xx, X3=xx, X4=xx)
  Y=list(Y1=matrix(scale(out_cohort1$n_krea), ncol=1), Y2=matrix(scale(out_cohort1$n_harn), ncol=1),  
         Y3=matrix(scale(out_cohort1$n_laktat), ncol=1), Y4=matrix(out_cohort1$event, ncol=1))
  Y[[4]][Y[[4]]==0,]=-1
  
  X2=dat_cohort2
  X2$c_dialyse=NULL
  X2$event=NULL
  X2=as.matrix(X2)
  X_validate = apply(X2, 2, function(x)scale(x))
  Y2 = dat_cohort2$event
  Y2[Y2==0]=-1
  
  opts <- list(init=0, tol=10^-6, maxIter=10000, ter=2)
  ctasks=4
  
  cvResult = CV_multi_outcomes(X=X, Y=Y, nfolds=10, nlambda=200, lam_ratio=0.001, C=5, C2=10, ctasks=4, opts=opts, intercept=F)
  fitMTLComb=MTLComb_Train(X=X, Y=Y, nlambda=200, lambda=cvResult$lambda.weighted.min, ctasks=4, opts=opts, C=5, C2=10, intercept=F)
  predScores=X_validate %*% fitMTLComb$ws[[1]][,4]
  predAucVec = c(predAucVec, as.numeric(pROC::auc(as.factor(Y2), predScores)))
  
  
  ################################################################################
  # Binarization approach 
  ################################################################################
  X1=dat_cohort1
  X1$c_dialyse=NULL
  X1$event=NULL
  X1=as.matrix(X1)
  xx=apply(X1, 2, function(x)scale(x))
  X=list(X1=xx, X2=xx, X3=xx, X4=xx)
  Y=list(Y1=matrix(sign(out_cohort1$n_krea-median(out_cohort1$n_krea)), ncol=1), Y2=matrix(sign(out_cohort1$n_harn-median(out_cohort1$n_harn)), ncol=1),  
         Y3=matrix(sign(out_cohort1$n_laktat-median(out_cohort1$n_laktat)), ncol=1), Y4=matrix(out_cohort1$event, ncol=1))
  Y=lapply(Y, function(x){x[x==0,]=-1; return(x)})
  
  X2=dat_cohort2
  X2$c_dialyse=NULL
  X2$event=NULL
  X2=as.matrix(X2)
  X_validate = apply(X2, 2, function(x)scale(x))
  Y2 = dat_cohort2$event
  Y2[Y2==0]=-1
  
  opts <- list(init=0, tol=10^-6, maxIter=10000, ter=2)
  
  #####################
  #test cross-validation
  #####################
  CVsplit = split( sample(c(1:nrow(X[[1]]))) , seq(1,nrow(X[[1]]),nrow(X[[1]])/10))
  auc_fold=vector()
  lam_seq = vector()
  ntasks=length(X)
  for(i in 1:length(CVsplit)){
    curtrain = unlist(CVsplit[-i])
    curtest  = unlist(CVsplit[i])
    Xtrain = lapply(X, function(x)x[curtrain,])
    Xtest = lapply(X, function(x)x[curtest,])
    Ytrain = lapply(Y, function(x)x[curtrain,, drop=F])
    Ytest = lapply(Y, function(x)x[curtest,, drop=F])
    fitTrain = MTL_L21_Train(X=Xtrain, Y=Ytrain, type="classify", nlambda=200, lam_ratio=0.001, C=10, opts=opts)
    
    yhatC=lapply(fitTrain$ws, function(w)lapply(1:ntasks, function(x) {score=exp(Xtest[[x]]%*%w[,x]); return(score/(1+score))}))
    auc_fold=rbind(auc_fold, sapply(yhatC, function(x)mean(sapply(1:ntasks, function(xx)as.numeric(pROC::auc(as.factor(Ytest[[xx]]), x[[xx]], quiet=TRUE))))))
    
    lam_seq = rbind(lam_seq, fitTrain$lam_seq)
  }
  lambda.min=colMeans(lam_seq)[order(colMeans(auc_fold), decreasing = T)[1]]
  
  fitMTLBin = MTL_L21_Train(X=X, Y=Y, type="classify", nlambda=200, lambda=lambda.min, C=10, opts=opts)
  predScoreMTLBin=X_validate %*% fitMTLBin$ws[[1]]
  predAucVec = c(predAucVec, as.numeric(pROC::auc(as.factor(Y2), predScoreMTLBin[,4])))
  
  coefMat = cbind(as.matrix(fitRidge$beta), as.matrix(fitLasso$beta), fitRF$importance, as.matrix(t(fitSVM$SV) %*% fitSVM$coefs),
                  as.matrix(t(fitSVM2$SV) %*% fitSVM2$coefs), as.matrix(t(fitSVM3$SV) %*% fitSVM3$coefs), 
                  fitMTLComb$ws[[1]][,4], fitMTLBin$ws[[1]][,4])
  return(list(predAucVec, coefMat))
}


saveRDS(file="1_results_rep_10_cohort_1_to_2_prediction.rds", rData)


