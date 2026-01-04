rm(list=ls())
gc()


#######################################################################################################################
#
#                             Cross-cohort prediction: Train on cohort 2 and test on cohort 1
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


predAucMat=vector()
coefList = list()
rData = foreach (rep = 1:10)%dopar%{
  print(rep)
  ################################################################################
  # Machine learning analysis on all variables
  ################################################################################
  X1=dat_cohort2
  X1$c_dialyse=NULL
  X1$event=NULL
  X1=as.matrix(X1)
  X1=apply(X1, 2, function(x)scale(x))
  Y1=dat_cohort2$event
  
  X2=dat_cohort1
  X2$c_dialyse=NULL
  X2$event=NULL
  X2=as.matrix(X2)
  X2 = apply(X2, 2, function(x)scale(x))
  Y2=dat_cohort1$event
  
  predAucVec = vector()
  coefMat = vector()
  
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
  X1=dat_cohort2
  X1$c_dialyse=NULL
  X1$event=NULL
  X1=as.matrix(X1)
  xx=apply(X1, 2, function(x)scale(x))
  X=list(X1=xx, X2=xx, X3=xx, X4=xx)
  Y=list(Y1=matrix(scale(out_cohort2$n_krea), ncol=1), Y2=matrix(scale(out_cohort2$n_harn), ncol=1),  
         Y3=matrix(scale(out_cohort2$n_laktat), ncol=1), Y4=matrix(out_cohort2$event, ncol=1))
  Y[[4]][Y[[4]]==0,]=-1
  
  X2=dat_cohort1
  X2$c_dialyse=NULL
  X2$event=NULL
  X2=as.matrix(X2)
  X_validate = apply(X2, 2, function(x)scale(x))
  Y2 = dat_cohort1$event
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
  X1=dat_cohort2
  X1$c_dialyse=NULL
  X1$event=NULL
  X1=as.matrix(X1)
  xx=apply(X1, 2, function(x)scale(x))
  X=list(X1=xx, X2=xx, X3=xx, X4=xx)
  Y=list(Y1=matrix(sign(out_cohort2$n_krea-median(out_cohort2$n_krea)), ncol=1), Y2=matrix(sign(out_cohort2$n_harn-median(out_cohort2$n_harn)), ncol=1),  
         Y3=matrix(sign(out_cohort2$n_laktat-median(out_cohort2$n_laktat)), ncol=1), Y4=matrix(out_cohort2$event, ncol=1))
  Y=lapply(Y, function(x){x[x==0,]=-1; return(x)})
  
  X2=dat_cohort1
  X2$c_dialyse=NULL
  X2$event=NULL
  X2=as.matrix(X2)
  X_validate = apply(X2, 2, function(x)scale(x))
  Y2 = dat_cohort1$event
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


saveRDS(file="2_results_rep_10_cohort_2_to_1_prediction.rds", rData)




rData1 = readRDS("1_results_rep_10_cohort_1_to_2_prediction.rds")
rData2 = readRDS("2_results_rep_10_cohort_2_to_1_prediction.rds")

rowMeans(sapply(rData1, function(x)x[[1]]))
#0.7325763 0.6809639 0.6733617 0.6530860 0.7233010 0.6867198 0.7021151 0.7002254

rowMeans(sapply(rData2, function(x)x[[1]]))
#0.7348265 0.7033537 0.7109131 0.6476704 0.7145091 0.6961382 0.7501563 0.7356942


rowMeans(sapply(1:10, function(x)sapply(1:8, function(y)cor(rData1[[x]][[2]][,y], rData2[[x]][[2]][,y]))))
#0.40571079 -0.02289282  0.65221077  0.21920604 -0.46533669 -0.54403207  0.69719999  0.58269743

numSharedMat=vector()
for (i in 1:10){
  coef1 = rData1[[i]][[2]]
  coef2 = rData2[[i]][[2]]
  sharedVec=vector()
  for (j in 1:8){
    v1 = coef1[,j]
    v2 = coef2[,j]
    num1=sum(v1!=0); 
    if (num1>10) {num1=10}
    num2=sum(v2!=0)
    if (num2>10) {num2=10}
    sharedVec=c(sharedVec, length(intersect(order(abs(v1), decreasing = T)[1:num1], order(abs(v2), decreasing = T)[1:num2])))
  }
  numSharedMat=rbind(numSharedMat, sharedVec)
}

Y2=list(Y1=matrix(scale(out_cohort2$n_krea), ncol=1), Y2=matrix(scale(out_cohort2$n_harn), ncol=1), Y3=matrix(scale(out_cohort2$n_laktat), ncol=1))
Y1=list(Y1=matrix(scale(out_cohort1$n_krea), ncol=1), Y2=matrix(scale(out_cohort1$n_harn), ncol=1), Y3=matrix(scale(out_cohort1$n_laktat), ncol=1))
X1=dat_cohort1
X1$c_dialyse=NULL
X1$event=NULL
X1=as.matrix(X1)
X2=dat_cohort2
X2$c_dialyse=NULL
X2$event=NULL
X2=as.matrix(X2)

coefList1=lapply(rData1, function(x)x[[2]])
predScore1 = lapply(coefList1, function(x)X2%*%x)
predScore1=Reduce("+", predScore1)/10
apply(predScore1, 2, function(x)sapply(Y2, function(y)summary(lm(y~x))$r.squared))
# Y1 0.01773935 0.002785872       0.02393678 0.012103142 0.02996903 0.03598388 0.04039412 0.02728710
# Y2 0.02986976 0.007459014       0.03159301 0.004008452 0.06722702 0.11640168 0.09329688 0.06365080
# Y3 0.07187036 0.041176339       0.01900051 0.018063749 0.10657342 0.08414494 0.12344276 0.09731231

coefList2=lapply(rData2, function(x)x[[2]])
predScore2 = lapply(coefList2, function(x)X1%*%x)
predScore2=Reduce("+", predScore2)/10
apply(predScore2, 2, function(x)sapply(Y1, function(y)summary(lm(y~x))$r.squared))
# Y1 0.007436608 0.015081569      0.004207778 0.0120021083 0.03670503 0.04650076 0.06768309 0.12122045
# Y2 0.047423023 0.053235485      0.020973330 0.0002408533 0.14375705 0.11747892 0.20682291 0.30593787
# Y3 0.028728119 0.009107434      0.007955716 0.0035351635 0.05320992 0.06451891 0.07414155 0.09156665

