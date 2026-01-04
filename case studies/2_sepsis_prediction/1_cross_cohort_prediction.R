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
plot(cvR)

fitRidge = glmnet(x=X1, y=Y1, family="binomial", lambda = cvR$lambda.min, alpha = 0)
fitRidge$beta[order(abs(fitRidge$beta), decreasing = T),]

predScore=predict(fitRidge, X2)
coef(summary(glm(as.factor(Y2)~predScore, family = binomial)))
#             Estimate Std. Error  z value     Pr(>|z|)
#predScore   1.2169025  0.2728765 4.459536 8.213722e-06
as.numeric(pROC::auc(as.factor(Y2), predScore))
#0.7319695

cvR=cv.glmnet(x=X1, y=Y1, nfolds = 10, family="binomial", alpha = 1)
fitLasso = glmnet(x=X1, y=Y1, family="binomial", lambda = cvR$lambda.min, alpha = 1)
predScoreLasso=predict(fitLasso, X2)
coef(summary(glm(as.factor(Y2)~predScoreLasso, family = binomial)))
#               Estimate Std. Error   z value     Pr(>|z|)
#predScoreLasso 0.90830120  0.2407894 3.7721813 0.0001618266
as.numeric(pROC::auc(as.factor(Y2), predScoreLasso))
#0.6820388

library(randomForest)
fitRF = randomForest(y = as.factor(Y1), x = X1, ntree=5000)
predScoreRF = predict(fitRF, X2, type="prob")[,2]
coef(summary(glm(as.factor(Y2)~predScoreRF, family = binomial)))
#           Estimate Std. Error   z value     Pr(>|z|)
#predScoreRF  3.926800  1.1616682  3.380311 7.240378e-04
as.numeric(pROC::auc(as.factor(Y2), predScoreRF))
#0.6725035

library(e1071)
fitSVM = svm(y=as.factor(Y1), x=X1, kernel = "linear", cost = 1, scale = FALSE, probability=TRUE)
predScoreSVM = attr(predict(fitSVM, X2,probability=TRUE), "probabilities")[,2]
coef(summary(glm(as.factor(Y2)~predScoreSVM, family = binomial)))
#           Estimate Std. Error   z value     Pr(>|z|)
#predScoreSVM  6.989357  2.1375782  3.269755 1.076407e-03
as.numeric(pROC::auc(as.factor(Y2), predScoreSVM))
#0.653086

fitSVM2 = svm(y=as.factor(Y1), x=X1, kernel = "radial", cost = 1, scale = FALSE, probability=TRUE)
predScoreSVM2 = attr(predict(fitSVM2, X2,probability=TRUE), "probabilities")[,2]
coef(summary(glm(as.factor(Y2)~predScoreSVM2, family = binomial)))
#               Estimate Std. Error   z value     Pr(>|z|)
#predScoreSVM2  4.717790  1.1267142  4.187211 2.824034e-05
as.numeric(pROC::auc(as.factor(Y2), predScoreSVM2))
#0.723301

fitSVM3 = svm(y=as.factor(Y1), x=X1, kernel = "polynomial", cost = 1, scale = FALSE, probability=TRUE)
predScoreSVM3 = attr(predict(fitSVM3, X2,probability=TRUE), "probabilities")[,2]
coef(summary(glm(as.factor(Y2)~predScoreSVM3, family = binomial)))
#               Estimate Std. Error   z value     Pr(>|z|)
#predScoreSVM3  4.577490   1.606121  2.850027 0.0043715460
as.numeric(pROC::auc(as.factor(Y2), predScoreSVM3))
#0.6867198
################################################################################








################################################################################
# MTLComb analysis with inductive CV
################################################################################
source('../../MTLComb/MTLComb.R')
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

#####################
#test cross-validation
#####################
cvResult = CV_multi_outcomes(X=X, Y=Y, nfolds=10, nlambda=200, lam_ratio=0.001, C2=10, C=5, ctasks=4, opts=opts, intercept=F)

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
#z value  3.282835057 3.5510786923 3.6158027520 4.018502e+00
#Pr(>|z|) 0.001027688 0.0003836557 0.0002994183 5.856922e-05
sapply(1:4, function(x)as.numeric(pROC::auc(as.factor(Y2), predScores[,x])))
#0.6690361 0.6792649 0.6835992 0.7102982

sort(sqrt(rowSums(fitMTLComb$ws[[1]]^2)), decreasing = T)
fitMTLComb$ws[[1]][names(sort(rowSums(fitMTLComb$ws[[1]]^2), decreasing = T)[1:4]),]
# n_sapsii 0.06133642 0.17208894 0.08133984 0.20338815
# n_ptt    0.07497106 0.18042299 0.09653156 0.16991894
# n_krea   0.18160193 0.07909906 0.04373691 0.02125854
# n_iss    0.08683392 0.07633916 0.08875936 0.12123054
################################################################################








################################################################################
# Binarization approach 
################################################################################
source('../../MTLComb/MTL_L21.R')
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
sapply(1:4, function(x)coef(summary(glm(as.factor(Y2)~predScoreMTLBin[,x], family=binomial)))[2,3:4])
#z value  3.2929078520 3.5139112607 3.4619719384 3.6918939158
#Pr(>|z|) 0.0009915697 0.0004415603 0.0005362331 0.0002225903

sapply(1:4, function(x)as.numeric(pROC::auc(as.factor(Y2), predScoreMTLBin[,x])))
#0.6622746 0.6777046 0.6744105 0.6848128
################################################################################


Models = list(fitRidge=fitRidge, fitMTLComb=fitMTLComb, fitMTLBin=fitMTLBin, fitLasso=fitLasso, fitRF=fitRF, fitSVM=fitSVM, fitSVM2=fitSVM2, fitSVM3=fitSVM3)
PredScores = list(predScore=predScore, predScores=predScores, predScoreMTLBin=predScoreMTLBin, predScoreLasso=predScoreLasso,
                  predScoreRF=predScoreRF, predScoreSVM=predScoreSVM, predScoreSVM2=predScoreSVM2, predScoreSVM3=predScoreSVM3)
cohort1_result=list(Models=Models, PredScores=PredScores)

saveRDS(file="1_cohort_1_to_2_data_result.rds", cohort1_result)





