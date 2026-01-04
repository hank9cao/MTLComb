rm(list=ls())
gc()


source ("../../../MTLComb/MTLComb.R")
library(foreach)
library(doParallel)
registerDoParallel(10)
library(randomForest)
library(pROC)
library(glmnet)


# data
x_Dx_1 = readRDS("../../preprocessed_data/HBCC/gene_expression_predict_Dx.rds")
y_Dx_1 = readRDS("../../preprocessed_data/HBCC/y_matrix_predict_Dx.rds")
x_age_1 = readRDS("../../preprocessed_data/HBCC/gene_expression_predict_age.rds")
y_age_1 = readRDS("../../preprocessed_data/HBCC/y_matrix_predict_age.rds")

x_Dx_2 = readRDS("../../preprocessed_data/GEO/gene_expression_predict_Dx.rds")
y_Dx_2 = readRDS("../../preprocessed_data/GEO/y_matrix_predict_Dx.rds")
x_age_2 = readRDS("../../preprocessed_data/GEO/gene_expression_predict_age.rds")
y_age_2 = readRDS("../../preprocessed_data/GEO/y_matrix_predict_age.rds")

# check if subjects are aligned 
dim(x_Dx_1) #359
dim(x_Dx_2) #317
sum(dim(rownames(x_Dx_1)!=rownames(y_Dx_1)))
sum(dim(rownames(x_Dx_2)!=rownames(y_Dx_2)))
sum(dim(rownames(x_age_1)!=rownames(y_age_1)))
sum(dim(rownames(x_age_2)!=rownames(y_age_2)))


#unify gene space
inteGenes = intersect(intersect(intersect(colnames(x_Dx_1), colnames(x_Dx_2)), colnames(x_age_1)), colnames(x_age_2)) #8799
x_Dx_1 = as.matrix(x_Dx_1[, match(inteGenes, colnames(x_Dx_1))])
x_Dx_2 = as.matrix(x_Dx_2[, match(inteGenes, colnames(x_Dx_2))])
x_age_1 = as.matrix(x_age_1[, match(inteGenes, colnames(x_age_1))])
x_age_2 = as.matrix(x_age_2[, match(inteGenes, colnames(x_age_2))])
sum(colnames(x_Dx_1)!=colnames(x_Dx_2))
sum(colnames(x_Dx_1)!=colnames(x_age_1))
sum(colnames(x_Dx_1)!=colnames(x_age_2))

# prepare outcomes
y_Dx_1 = matrix(2*(as.numeric(y_Dx_1[,"diagnosis"])-0.5), ncol=1)
rownames(y_Dx_1)=rownames(x_Dx_1)
y_Dx_2 = 2*(matrix(as.numeric(y_Dx_2[,"label"]), ncol=1)-1.5)
rownames(y_Dx_2)=rownames(x_Dx_2)
y_age_2=matrix(y_age_2[,"age"], ncol=1)
rownames(y_age_2)=rownames(x_age_2)
y_age_1=matrix(y_age_1[,"age"], ncol=1)
rownames(y_age_1)=rownames(x_age_1)


# Prepare data lists for analysis
X1=list(scale(x_Dx_1), scale(x_age_1))
X2=list(scale(x_Dx_2), scale(x_age_2))
Y1=list(y_Dx_1, scale(y_age_1))
Y2=list(y_Dx_2, scale(y_age_2))

#hyperparameters
ctasks=1
opts=list(init=0, maxIter=200, tol=0.001, ter=2)




#######################################################################################
# Train on cohort 1 and test on cohort 2
#######################################################################################
#MTLComb
cvFit = MTLComb_CV(X=X1, Y=Y1, nfolds=10, nlambda=10, lam_ratio=0.01, ctasks=ctasks, opts=opts, C=5, C2=1, intercept=F)
fit=MTLComb_Train(X=X1, Y=Y1, nlambda=20, lambda=cvFit$lambda.weighted.min, ctasks=ctasks, opts=opts, C=5, C2=1, intercept=F)
plot(fit$ws[[1]][rowSums(fit$ws[[1]])!=0,])
dim(fit$ws[[1]][rowSums(fit$ws[[1]])!=0,])
predScores = sapply(X2, function(x)x%*%rowMeans(fit$ws[[1]]))
coef(summary(glm(as.factor(Y2[[1]])~predScores[[1]], family = binomial)))
#predScores[[1]] 12.84762363  3.1429077 4.0878145 4.354562e-05
roc(as.factor(Y2[[1]]), predScores[[1]])
#AUC=0.7077
summary(lm(Y2[[2]]~predScores[[2]]))
#predScores[[2]] 5.861e+00  1.913e+00   3.063  0.00289 **
#R2: 0.08435   

saveRDS(file="fit_test.rds", fit)

#build model
X=list(X1[[1]], X2[[1]], X1[[2]], X2[[2]])
Y=list(Y1[[1]], Y2[[1]], Y1[[2]], Y2[[2]])
ctasks=c(1,2)
cvFit_all = MTLComb_CV(X=X, Y=Y, nfolds=10, nlambda=10, lam_ratio=0.01, ctasks=ctasks, opts=opts, C=10, C2=10, intercept=F)
fit_all = MTLComb_Train(X=X, Y=Y, nlambda=20, lambda=cvFit_all$lambda.weighted.min, ctasks=ctasks, opts=opts, C=10, C2=10, intercept=F)

saveRDS(file="fit_cvFit_all.rds", list(fit=fit_all, cvFit=cvFit_all))







#################################################################################################################################
#Analyze models
fits=readRDS("fit1_fit2_test.rds")
fit1=fits[[1]]
fit2=fits[[2]]
selectedGenes1 = names(which(rowSums(fit1$ws[[1]])!=0))
selectedGenes2 = names(which(rowSums(fit2$ws[[1]])!=0))
num1 = length(selectedGenes1)
num2 = length(selectedGenes2)
topgenes1 = names(sort(abs(rowSums(fit1$ws[[1]])[selectedGenes1]), decreasing = T)[1:(0.25*num1)])
topgenes2 = names(sort(abs(rowSums(fit2$ws[[1]])[selectedGenes2]), decreasing = T)[1:(0.25*num2)])
intersect(topgenes1, topgenes2)
#[1] "GALNT16" "KCNK1"   "DACH2"   "ANGPTL2" "CACNA1G" "PYGB"    "RHBDL3"  "ADO"     "IGDCC4"  "DDX54"   "LPL"     "CA11"    "DUSP3"  
#[14] "PRAG1"   "FAM81A"  "PHF23"   "HOMER1"  "SBK1"    "NETO2"   "TGFB3"  
topgenes1 = names(sort(abs(rowSums(fit1$ws[[1]])[selectedGenes1]), decreasing = T)[1:100])
topgenes2 = names(sort(abs(rowSums(fit2$ws[[1]])[selectedGenes2]), decreasing = T)[1:100])
intersect(topgenes1, topgenes2)
#[1] "GALNT16"  "IGDCC4"   "LPL"      "PRAG1"    "C1orf52"  "ZMPSTE24"





#################
#RF
#################
fitRF = randomForest(y = as.factor(Y1[[1]]), x = X1[[1]], ntree=5000)
saveRDS(file="RFmodel_Dx_cohort1", fitRF)
fitRF=readRDS("RFmodel_Dx_cohort1")
predScoreRF = predict(fitRF, X2[[1]], type="prob")[,2]
coef(summary(glm(as.factor(Y2[[1]])~predScoreRF, family = binomial)))
#predScoreRF 15.217863  1.9010521  8.004969 1.194966e-15
as.numeric(pROC::auc(as.factor(Y2[[1]]), predScoreRF))
#AUC=0.7934147
predScoreRF = predict(fitRF, X2[[2]], type="prob")[,2]
summary(lm(Y2[[2]]~predScoreRF))
#predScoreRF  -2.9357     0.5943  -4.940 1.34e-06 ***
#R2=0.07588 

fitRF = randomForest(y = Y1[[2]], x = X1[[2]], ntree=5000)
saveRDS(file="RFmodel_Age_cohort1", fitRF)
fitRF=readRDS("RFmodel_Age_cohort1")
predScoreRF = predict(fitRF, X2[[1]])
coef(summary(glm(as.factor(Y2[[1]])~predScoreRF, family = binomial)))
#predScoreRF -0.7014419  0.6627273 -1.058417 0.2898653643
as.numeric(pROC::auc(as.factor(Y2[[1]]), predScoreRF))
#AUC=0.5380662
predScoreRF = predict(fitRF, X2[[2]])
summary(lm(Y2[[2]]~predScoreRF))
#predScoreRF  2.11194    0.15051  14.032   <2e-16 ***
#R2: 0.4073 

fitRF = randomForest(y = as.factor(Y2[[1]]), x = X2[[1]], ntree=5000)
saveRDS(file="RFmodel_Dx_cohort2", fitRF)
fitRF=readRDS("RFmodel_Dx_cohort2")
predScoreRF = predict(fitRF, X1[[1]], type="prob")[,2]
coef(summary(glm(as.factor(Y1[[1]])~predScoreRF, family = binomial)))
#predScoreRF  9.232726  1.2404919  7.442794 9.857741e-14
as.numeric(pROC::auc(as.factor(Y1[[1]]), predScoreRF))
#AUC=0.7466465
predScoreRF = predict(fitRF, X1[[2]])
summary(lm(Y1[[2]]~predScoreRF))
#predScoreRF1  0.06564    0.15863   0.414    0.679
#R2=-0.002948 

fitRF = randomForest(y = Y2[[2]], x = X2[[2]], ntree=5000)
saveRDS(file="RFmodel_Age_cohort2", fitRF)
fitRF=readRDS("RFmodel_Age_cohort2")
predScoreRF = predict(fitRF, X1[[1]])
coef(summary(glm(as.factor(Y1[[1]])~predScoreRF, family = binomial)))
#predScoreRF -1.1327796  0.5122276 -2.211477 0.027002823
as.numeric(pROC::auc(as.factor(Y1[[1]]), predScoreRF))
#AUC=0.5592921
predScoreRF = predict(fitRF, X1[[2]])
summary(lm(Y1[[2]]~predScoreRF))
#predScoreRF  1.60972    0.21435   7.510 7.90e-13 ***
#R2=0.1642 


#################
#Ridge
#################
cvR=cv.glmnet(x=X1[[1]], y=as.factor(Y1[[1]]), nfolds = 10, family="binomial", alpha = 0)
fitRidge = glmnet(x=X1[[1]], y=as.factor(Y1[[1]]), family="binomial", lambda = cvR$lambda.min, alpha = 0)
predScore=predict(fitRidge, X2[[1]])
pROC::auc(as.factor(Y2[[1]]), predScore[,1])
#0.7512
predScore=predict(fitRidge, X2[[2]])
summary(lm(Y2[[2]]~predScore))
#predScore     0.9786     0.3596   2.721  0.00691 **
#r2=0.02198 

cvR=cv.glmnet(x=X1[[2]], y=Y1[[2]], nfolds = 10, alpha = 0)
fitRidge = glmnet(x=X1[[2]], y=Y1[[2]], lambda = cvR$lambda.min, alpha = 0)
predScore=predict(fitRidge, X2[[1]])
pROC::auc(as.factor(Y2[[1]]), predScore[,1])
#auc=0.5781
predScore=predict(fitRidge, X2[[2]])
summary(lm(Y2[[2]]~predScore))
#predScore   1.919e+00  3.552e-01   5.402  1.4e-07 ***
#R2=0.08998 

cvR=cv.glmnet(x=X2[[1]], y=as.factor(Y2[[1]]), nfolds = 10, family="binomial", alpha = 0)
fitRidge = glmnet(x=X2[[1]], y=as.factor(Y2[[1]]), family="binomial", lambda = cvR$lambda.min, alpha = 0)
predScore=predict(fitRidge, X1[[1]])
pROC::auc(as.factor(Y1[[1]]), predScore[,1])
#0.7306
predScore=predict(fitRidge, X1[[2]])
summary(lm(Y1[[2]]~predScore))
#predScore     0.6013     0.2329   2.582   0.0103 *
#0.01969 

cvR=cv.glmnet(x=X2[[2]], y=Y2[[2]], nfolds = 10, alpha = 0)
fitRidge = glmnet(x=X2[[2]], y=Y2[[2]], lambda = cvR$lambda.min, alpha = 0)
predScore=predict(fitRidge, X1[[1]])
pROC::auc(as.factor(Y1[[1]]), predScore[,1])
#0.5889
predScore=predict(fitRidge, X1[[2]])
coef(summary(lm(Y1[[2]]~predScore)))
#predScore    4.142128e+00 0.53454372  7.748903e+00 1.700449e-13
#R2=0.1731


#################
#Lasso
#################
cvR=cv.glmnet(x=X1[[1]], y=as.factor(Y1[[1]]), nfolds = 10, family="binomial", alpha = 1)
fitLasso = glmnet(x=X1[[1]], y=as.factor(Y1[[1]]), family="binomial", lambda = cvR$lambda.min, alpha = 1)
predScore=predict(fitLasso, X2[[1]])
pROC::auc(as.factor(Y2[[1]]), predScore[,1])
#0.6854
cvR=cv.glmnet(x=X1[[2]], y=Y1[[2]], nfolds = 10, alpha = 1)
fitLasso = glmnet(x=X1[[2]], y=Y1[[2]], lambda = cvR$lambda.min, alpha = 1)
predScore=predict(fitLasso, X2[[2]])
coef(summary(lm(Y2[[2]]~predScore)))
#predScore    3.044975e+00 0.75560293  4.029861e+00 7.172492e-05

cvR=cv.glmnet(x=X2[[1]], y=as.factor(Y2[[1]]), nfolds = 10, family="binomial", alpha = 1)
fitLasso = glmnet(x=X2[[1]], y=as.factor(Y2[[1]]), family="binomial", lambda = cvR$lambda.min, alpha = 1)
sum(fitLasso$beta!=0) #33
predScore=predict(fitLasso, X1[[1]])
pROC::auc(as.factor(Y1[[1]]), predScore[,1])
#0.6903
cvR=cv.glmnet(x=X2[[2]], y=Y2[[2]], nfolds = 10, alpha = 1)
fitLasso = glmnet(x=X2[[2]], y=Y2[[2]], lambda = cvR$lambda.min, alpha = 1)
sum(fitLasso$beta!=0) #0




