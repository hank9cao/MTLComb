rm(list=ls())
gc()


source ("../../../MTLComb/MTLComb.R")
library(foreach)
library(doParallel)
registerDoParallel(10)
library(pROC)

# data
x_data1 = readRDS("../../preprocessed_data/HBCC/gene_expression_predict_Dx.rds")
y_data1 = readRDS("../../preprocessed_data/HBCC/y_matrix_predict_Dx.rds")
x_data2 = readRDS("../../preprocessed_data/HBCC/gene_expression_predict_age.rds")
y_data2 = readRDS("../../preprocessed_data/HBCC/y_matrix_predict_age.rds")


# harmonize the data 
head(y_data1)
x_data1[1:10,1:5]
sum(dim(rownames(y_data1)!=rownames(x_data1)))
dim(x_data1) #359
X1 = as.matrix(x_data1)
Y1 = matrix(2*(as.numeric(y_data1[,"diagnosis"])-0.5), ncol=1)
rownames(Y1)=rownames(y_data1)
table(Y1)

head(y_data2)
x_data2[1:10,1:5]
sum(dim(rownames(y_data2)!=rownames(x_data2)))
dim(x_data2) #283
X2 = as.matrix(x_data2)
Y2 = matrix(y_data2[,"age"], ncol=1)
rownames(Y2)=rownames(y_data2)
str(Y2)


# Prepare data lists for nested CV
inte = intersect(colnames(X1), colnames(X2)) #9569
X1 = X1[, match(inte, colnames(X1))]
X2 = X2[, match(inte, colnames(X2))]
X=list(as.matrix(scale(X1)), as.matrix(scale(X2)))
Y=list(Y1, scale(Y2))


# Nested CV start
set.seed(1123213)
CVsplit = lapply(c(1,2), function(x)split( sample(c(1:nrow(X[[x]]))) , seq(1,nrow(X[[x]]),nrow(X[[x]])/10)))
r2_Ave=0
lam_seq = vector()
ntasks=length(X)
rtasks=2
ctasks=1
opts=list(init=0, maxIter=200, tol=0.001, ter=2)
rData = foreach(i = 1:10, .combine=c)%dopar%{
  print(i)
  
  Xtrain=list(); Xtest=list(); Ytrain=list(); Ytest=list();
  for (j in 1:length(CVsplit)){
    curtrain = unlist(CVsplit[[j]][-i])
    curtest  = unlist(CVsplit[[j]][i])
    Xtrain[[j]] = X[[j]][curtrain,]
    Xtest[[j]] = X[[j]][curtest,]
    Ytrain[[j]] = Y[[j]][curtrain,, drop=F]
    Ytest[[j]] = Y[[j]][curtest,, drop=F]
  }
  
  cvFit = MTLComb_CV(X=Xtrain, Y=Ytrain, nfolds=10, nlambda=10, lam_ratio=0.01, ctasks=ctasks, opts=opts, C=5, C2=1, intercept=F)
  
  fit=MTLComb_Train(X=Xtrain, Y=Ytrain, nlambda=10, lambda=cvFit$lambda.weighted.min, ctasks=ctasks, opts=opts, C=5, C2=1, intercept=F)
  
  predScores1 = Xtest[[1]]%*%rowMeans(fit$ws[[1]])
  predScores2 = Xtest[[2]]%*%rowMeans(fit$ws[[1]])
  
  return(list(cbind(predScores1, Ytest[[1]]), cbind(predScores2, Ytest[[2]])))
}
saveRDS(file="results_nested_cv.rds", rData)


rData_dx = do.call(rbind, sapply(seq(1,20,2), function(x)rData[[x]]))
rData_age = do.call(rbind, sapply(seq(2,20,2), function(x)rData[[x]]))

summary(lm(rData_age[,2]~rData_age[,1]))
#rData_age[, 1]  1.109810   0.391743   2.833  0.00495 **
#R2=0.02431 
coef(summary(glm(as.factor(rData_dx[,2])~rData_dx[,1], family = binomial)))
#rData_dx[, 1]   2.7468     0.8791   3.125  0.00178 **
roc(as.factor(rData_dx[,2]), rData_dx[,1])
#AUC=0.6377

cvFit = MTLComb_CV(X=X, Y=Y, nfolds=10, nlambda=10, lam_ratio=0.01, ctasks=ctasks, opts=opts, C=1, C2=1, intercept=F)
fit=MTLComb_Train(X=X, Y=Y, nlambda=10, lambda=cvFit$lambda.weighted.min, ctasks=ctasks, opts=opts, C=1, C2=1, intercept=F)
plot(fit$ws[[1]][,1], fit$ws[[1]][,2])
sort(rowSums(fit$ws[[1]][abs(rowSums(fit$ws[[1]]))!=0, ]))
saveRDS(file="fit_C=1_C2=1.rds", fit)











