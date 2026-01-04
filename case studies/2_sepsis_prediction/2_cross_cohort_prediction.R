rm(list=ls())
gc()


#######################################################################################################################
#
#                             Cross-cohort prediction: Train on cohort 1 and test on cohort 2
#
#######################################################################################################################



dat_cohort1 = readRDS("0_data_cohort1.rds")
dat_cohort2 = readRDS("0_data_cohort2.rds")

out_cohort1 = readRDS("0_out_cohort1.rds")
out_cohort2 = readRDS("0_out_cohort2.rds")


################################################################################
# Machine learning analysis on all variables
################################################################################
library(glmnet)
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

cvR=cv.glmnet(x=X1, y=Y1, nfolds = 10, family="binomial", alpha = 0)
plot(cvR)

fitRidge = glmnet(x=X1, y=Y1, family="binomial", lambda = cvR$lambda.min, alpha = 0)
fitRidge$beta[order(abs(fitRidge$beta), decreasing = T),]

predScore=predict(fitRidge, X2)
coef(summary(glm(as.factor(Y2)~predScore, family = binomial)))
#predScore    0.8104070  0.1524768  5.314954 1.066844e-07
as.numeric(pROC::auc(as.factor(Y2), predScore))
#0.7280331

cvR=cv.glmnet(x=X1, y=Y1, nfolds = 10, family="binomial", alpha = 1)
fitLasso = glmnet(x=X1, y=Y1, family="binomial", lambda = cvR$lambda.min, alpha = 1)
predScoreLasso=predict(fitLasso, X2)
coef(summary(glm(as.factor(Y2)~predScoreLasso, family = binomial)))
#predScoreLasso  0.8796384  0.1756033  5.009235 5.464680e-07
as.numeric(pROC::auc(as.factor(Y2), predScoreLasso))
#0.7035647

library(randomForest)
fitRF = randomForest(y = as.factor(Y1), x = X1, ntree=5000)
predScoreRF = predict(fitRF, X2, type="prob")[,2]
coef(summary(glm(as.factor(Y2)~predScoreRF, family = binomial)))
#predScoreRF  5.414409  1.0782191  5.021622 5.123692e-07
as.numeric(pROC::auc(as.factor(Y2), predScoreRF))
#0.7131801

library(e1071)
fitSVM = svm(y=as.factor(Y1), x=X1, kernel = "linear", cost = 1, scale = FALSE, probability=TRUE)
predScoreSVM = attr(predict(fitSVM, X2,probability=TRUE), "probabilities")[,2]
coef(summary(glm(as.factor(Y2)~predScoreSVM, family = binomial)))
#predScoreSVM -4.235293  1.1138073 -3.802537 0.0001432219
as.numeric(pROC::auc(as.factor(Y2), predScoreSVM))
#0.6476704

fitSVM2 = svm(y=as.factor(Y1), x=X1, kernel = "radial", cost = 1, scale = FALSE, probability=TRUE)
predScoreSVM2 = attr(predict(fitSVM2, X2,probability=TRUE), "probabilities")[,2]
coef(summary(glm(as.factor(Y2)~predScoreSVM2, family = binomial)))
#predScoreSVM2 -4.095095  0.8028592 -5.100639 3.385092e-07
as.numeric(pROC::auc(as.factor(Y2), predScoreSVM2))
#0.7145091

fitSVM3 = svm(y=as.factor(Y1), x=X1, kernel = "polynomial", cost = 1, scale = FALSE, probability=TRUE)
predScoreSVM3 = attr(predict(fitSVM3, X2,probability=TRUE), "probabilities")[,2]
coef(summary(glm(as.factor(Y2)~predScoreSVM3, family = binomial)))
#predScoreSVM3 -6.212980  1.5175328 -4.094132 4.237523e-05
as.numeric(pROC::auc(as.factor(Y2), predScoreSVM3))
#0.6961382
################################################################################








################################################################################
# MTLComb analysis with inductive CV
################################################################################
source('../../MTLComb/MTLComb.R')
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

#####################
#test cross-validation
#####################
cvResult = CV_multi_outcomes(X=X, Y=Y, nfolds=10, nlambda=200, lam_ratio=0.001, C=5, C2=10, ctasks=4, opts=opts, intercept=F)

####################
#training
####################
fitMTLComb=MTLComb_Train(X=X, Y=Y, nlambda=200, lambda=cvResult$lambda.weighted.min, ctasks=4, opts=opts, C=5, C2=10, intercept=F)
plot(sqrt(rowSums((fitMTLComb$ws[[1]])^2)))

#####################
#test prediction
#####################
predScores=X_validate %*% fitMTLComb$ws[[1]]
pairs(predScores)
sapply(1:4, function(x)coef(summary(glm(as.factor(Y2)~predScores[,x], family=binomial)))[2,3:4])
#z value  4.803498e+00 4.981445e+00 5.355252e+00 5.740914e+00
#Pr(>|z|) 1.559173e-06 6.311130e-07 8.543740e-08 9.416694e-09
sapply(1:4, function(x)as.numeric(pROC::auc(as.factor(Y2), predScores[,x])))
#0.7042683 0.7077861 0.7288931 0.7541432

sort(sqrt(rowSums(fitMTLComb$ws[[1]]^2)), decreasing = T)
fitMTLComb$ws[[1]][names(sort(rowSums(fitMTLComb$ws[[1]]^2), decreasing = T)[1:4]),]
# n_krea       0.03818995 0.02700299 0.010113030 0.012685814
# n_harn       0.02831975 0.03419138 0.016736303 0.007927592
# o_sofa_renal 0.03577605 0.02495304 0.009957463 0.007652806
# n_sapsii     0.01562637 0.02836763 0.016235231 0.015068901
################################################################################








################################################################################
# Binarization approach 
################################################################################
source('../../MTLComb/MTL_L21.R')
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
sapply(1:4, function(x)coef(summary(glm(as.factor(Y2)~predScoreMTLBin[,x], family=binomial)))[2,3:4])
#z value  4.391884e+00 4.797313e+00 4.722487e+00 5.655388e+00
#Pr(>|z|) 1.123726e-05 1.608082e-06 2.329783e-06 1.554949e-08

sapply(1:4, function(x)as.numeric(pROC::auc(as.factor(Y2), predScoreMTLBin[,x])))
#0.6774547 0.6926986 0.6874609 0.7453877
################################################################################


Models = list(fitRidge=fitRidge, fitMTLComb=fitMTLComb, fitMTLBin=fitMTLBin, fitLasso=fitLasso, fitRF=fitRF, fitSVM=fitSVM, fitSVM2=fitSVM2, fitSVM3=fitSVM3)
PredScores = list(predScore=predScore, predScores=predScores, predScoreMTLBin=predScoreMTLBin, predScoreLasso=predScoreLasso,
                  predScoreRF=predScoreRF, predScoreSVM=predScoreSVM, predScoreSVM2=predScoreSVM2, predScoreSVM3=predScoreSVM3)
cohort1_result=list(Models=Models, PredScores=PredScores)

saveRDS(file="2_cohort_2_to_1_data_result.rds", cohort1_result)





