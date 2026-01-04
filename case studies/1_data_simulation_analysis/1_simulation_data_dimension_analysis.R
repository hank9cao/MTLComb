rm(list=ls())
gc()

source('./MTLComb.R')
source('./MTL_L21.R')
library(foreach)
library(doParallel)
library(fmsb)
library(glmnet)
registerDoParallel(10)
library(glmnet)
library(randomForest)
library(e1071)

########################
# Main Procedure
#0, define data simulation functions
#1, run analysis
#2, collect resuslts
#3, plot results
########################







########################
#0, define data simulation functions
########################
#t=10; p=200; n=100; ctasks=1:5; noise=0.1; props=0.1
simulateData4MTLComb = function(t=10, p=200, n=100, props=41, ctasks=1:5, noise=0.1){
  rtasks=setdiff(1:t, ctasks)
  W <- matrix(data=sign(rnorm(t*p))*rnorm(t*p,mean=1),ncol=t, nrow = p)
  W[props:p,] <- 0
  X <- list(); Y <- list(); tX <- list(); tY <- list()
  for(i in 1:t){
    if (sum(i==ctasks)==1){
      X[[i]] <- matrix(data=rnorm(n*p),ncol=p, nrow = n)
      Y[[i]] <- sign(X[[i]] %*% W[,i] + noise * rnorm(n))
      tX[[i]] <- matrix(rnorm(p*n),nrow=n)
      tY[[i]] <- sign(tX[[i]] %*% W[, i] + rnorm(n) * noise)
    }else {
      X[[i]] <- matrix(data=rnorm(n*p),ncol=p, nrow = n)
      Y[[i]] <- X[[i]] %*% W[,i] + noise * rnorm(n)
      tX[[i]] <- matrix(rnorm(p*n),nrow=n)
      tY[[i]] <- tX[[i]] %*% W[, i] + rnorm(n) * noise
    }
  }
  X=lapply(X, function(x)scale(x))
  tX=lapply(tX, function(x)scale(x))
  Y[rtasks]=lapply(Y[rtasks], function(x)scale(x))
  Y[rtasks]=lapply(Y[rtasks], function(x)scale(x))
  simuData = list(X=X, Y=Y, tX=tX, tY=tY, W=W, ctasks=ctasks)
  return(simuData)
}



opts <- list(init=0, tol=10^-3, maxIter=500, ter=2)
propSubs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

########################
#1, run analysis
########################
resultDataRepeats = foreach(rep = 1:10)%dopar%{
  resultsData = list()
  for (i in 1:length(propSubs)){
    print(i)
    
    rData = simulateData4MTLComb(t=20, p=500, n=propSubs[i]*500, props=51, ctasks=1:10, noise=0.5)
    
    nSubs = dim(rData$X[[1]])[1]
    X = rData$X
    Y = rData$Y
    tX = rData$tX
    tY = rData$tY
    ntasks=length(X)
    ctasks=rData$ctasks
    rtasks=setdiff(1:ntasks, ctasks)
    numSigFeat = sum(rData$W[,1]!=0)
    sigFeatIdx = which(rData$W[,1]!=0)
    
    #############
    #MTLComb
    #############
    cvResult = MTLComb_CV(X=X, Y=Y, nfolds=5, lam_ratio=0.01, nlambda=200, ctasks=ctasks, opts=opts, C2=0.1, stratify = T)
    fit1 = MTLComb_Train(X=X, Y=Y, nlambda=200, lambda=cvResult$lambda.weighted.min, ctasks=ctasks, opts=opts, C2=0.1)
    predScores = sapply(1:ntasks, function(x)tX[[x]]%*%fit1$ws[[1]][,x])
    expVar1 = sapply(1:ntasks, function(x){
      if (sum(x == ctasks)==1){
        ev = NagelkerkeR2(glm(as.factor(tY[[x]])~predScores[,x], family = binomial))$R2
      } else {
        ev = summary(lm(tY[[x]]~predScores[,x]))$r.squared
      }
      return(ev)
    })
    numf=numSigFeat
    if (sum(sqrt(rowSums((fit1$ws[[1]])^2))!=0)<numSigFeat) {numf=sum(sqrt(rowSums((fit1$ws[[1]])^2))!=0)}
    featSel1 = length(intersect(order(sqrt(rowSums((fit1$ws[[1]])^2)), decreasing = T)[1:numf], sigFeatIdx))/numSigFeat
    
    
    
    #############
    #MTLBin
    #############
    bY = Y
    bY[rtasks] = lapply(bY[rtasks], function(x)matrix(sign(x-median(x))))
    bY=lapply(bY, function(x){x[x==0,]=-1; return(x)})
    cvResult = MTL_L21_CVInSite(X=X, Y=bY, type="classify", nfolds=5, lam_ratio=0.01, nlambda=200, opts=opts, stratify = T)
    fit2 = MTL_L21_Train(X=X, Y=bY, type="classify", nlambda=200, lambda=cvResult$lambda.min, opts=opts)
    predScores = sapply(1:ntasks, function(x)tX[[x]]%*%fit2$ws[[1]][,x])
    expVar2 = sapply(1:ntasks, function(x){
      if (sum(x == ctasks)==1){
        ev = NagelkerkeR2(glm(as.factor(tY[[x]])~predScores[,x], family = binomial))$R2
      } else {
        ev = summary(lm(tY[[x]]~predScores[,x]))$r.squared
      }
      return(ev)
    })
    numf=numSigFeat
    if (sum(sqrt(rowSums((fit2$ws[[1]])^2))!=0)<numSigFeat) {numf=sum(sqrt(rowSums((fit2$ws[[1]])^2))!=0)}
    featSel2 = length(intersect(order(sqrt(rowSums((fit2$ws[[1]])^2)), decreasing = T)[1:numf], sigFeatIdx))/numSigFeat
    
    
    #############
    #ridge, lasso, rf, svm
    #############
    predScoreMat1=vector()
    predScoreMat2=vector()
    predScoreMat3=vector()
    predScoreMat4=vector()
    predScoreMat5=vector()
    predScoreMat6=vector()
    betas = vector()
    for (j in 1:ntasks){
      if (sum(j == ctasks)==1){
        fam = "binomial"
        prob=TRUE
        yy=as.factor(tY[[j]])
      } else {
        fam = "gaussian"
        yy=tY[[j]]
        prob=FALSE
      }
      
      cvResult = cv.glmnet(x=X[[j]], y=yy, alpha=0, nfolds=5, family=fam )
      fit3 = glmnet(x=X[[j]], y=Y[[j]], alpha=0, lambda=cvResult$lambda.min, family=fam )
      predScoreMat1 = cbind(predScoreMat1, predict(fit3, tX[[j]]))
        
      cvResult = cv.glmnet(x=X[[j]], y=yy, alpha=1, nfolds=5, family=fam )
      fit4 = glmnet(x=X[[j]], y=Y[[j]], alpha=1, lambda=cvResult$lambda.min, family=fam )
      predScoreMat2 = cbind(predScoreMat2, predict(fit3, tX[[j]]))
      
      fit5 = randomForest(y = yy, x = X[[j]], ntree=5000)

      fit6 = svm(y=yy, x=X[[j]], kernel = "linear", cost = 1, scale = FALSE, probability=prob)

      fit7 = svm(y=yy, x=X[[j]], kernel = "radial", cost = 1, scale = FALSE, probability=prob)

      fit8 = svm(y=yy, x=X[[j]], kernel = "polynomial", cost = 1, scale = FALSE, probability=prob)
      
      if (sum(j == ctasks)==1){
        predScoreMat3 = cbind(predScoreMat3, predict(fit5, tX[[j]], type="prob")[,2])
        predScoreMat4 = cbind(predScoreMat4, attr(predict(fit6, tX[[j]], probability=prob), "probabilities")[,2])
        predScoreMat5 = cbind(predScoreMat5, attr(predict(fit7, tX[[j]], probability=prob), "probabilities")[,2])
        predScoreMat6 = cbind(predScoreMat6, attr(predict(fit8, tX[[j]], probability=prob), "probabilities")[,2])
      } else {
        predScoreMat3 = cbind(predScoreMat3, predict(fit5, tX[[j]]))
        predScoreMat4 = cbind(predScoreMat4, predict(fit6, tX[[j]], probability=prob))
        predScoreMat5 = cbind(predScoreMat5, predict(fit7, tX[[j]], probability=prob))
        predScoreMat6 = cbind(predScoreMat6, predict(fit8, tX[[j]], probability=prob))
      }
    }   
    
    expVar3 = sapply(ctasks, function(x){
      ev=vector()
      ev = c(ev, NagelkerkeR2(glm(as.factor(tY[[x]])~predScoreMat1[,x], family = binomial))$R2)
      ev = c(ev, NagelkerkeR2(glm(as.factor(tY[[x]])~predScoreMat2[,x], family = binomial))$R2)
      ev = c(ev, NagelkerkeR2(glm(as.factor(tY[[x]])~predScoreMat3[,x], family = binomial))$R2)
      ev = c(ev, NagelkerkeR2(glm(as.factor(tY[[x]])~predScoreMat4[,x], family = binomial))$R2)
      ev = c(ev, NagelkerkeR2(glm(as.factor(tY[[x]])~predScoreMat5[,x], family = binomial))$R2)
      ev = c(ev, NagelkerkeR2(glm(as.factor(tY[[x]])~predScoreMat6[,x], family = binomial))$R2)
      return(ev)
    })
    expVar3 = cbind(expVar3, sapply(rtasks, function(x){
      ev=vector()
      ev = c(ev, summary(lm(tY[[x]]~predScoreMat1[,x]))$r.squared)
      ev = c(ev, summary(lm(tY[[x]]~predScoreMat2[,x]))$r.squared)
      ev = c(ev, summary(lm(tY[[x]]~predScoreMat3[,x]))$r.squared)
      ev = c(ev, summary(lm(tY[[x]]~predScoreMat4[,x]))$r.squared)
      ev = c(ev, summary(lm(tY[[x]]~predScoreMat5[,x]))$r.squared)
      ev = c(ev, summary(lm(tY[[x]]~predScoreMat6[,x]))$r.squared)
      return(ev)
    }))
  
    featSel3 = vector()
    featSel3 = c(featSel3, length(intersect(order(abs(fit3$beta), decreasing = T)[1:numSigFeat], sigFeatIdx))/numSigFeat)
    
    numf = numSigFeat
    if (sum(abs(fit4$beta)!=0)<numSigFeat) {numf=sum(abs(fit4$beta)!=0)}
    if (numf==0){featSel3 = c(featSel3, 0)
    } else {featSel3 = c(featSel3, length(intersect(order(abs(fit4$beta), decreasing = T)[1:numf], sigFeatIdx))/numSigFeat)}
    
    featSel3 = c(featSel3, length(intersect(order(as.matrix(fit5$importance), decreasing = T)[1:numSigFeat], sigFeatIdx))/numSigFeat)
    
    featSel3 = c(featSel3, length(intersect(order(abs(as.matrix(t(fit6$SV) %*% fit6$coefs)), decreasing = T)[1:numSigFeat], sigFeatIdx))/numSigFeat)
    
    featSel3 = c(featSel3, length(intersect(order(abs(as.matrix(t(fit7$SV) %*% fit7$coefs)), decreasing = T)[1:numSigFeat], sigFeatIdx))/numSigFeat)
    
    featSel3 = c(featSel3, length(intersect(order(abs(as.matrix(t(fit8$SV) %*% fit8$coefs)), decreasing = T)[1:numSigFeat], sigFeatIdx))/numSigFeat)
    
    expVar = rbind(expVar1, expVar2, expVar3)
    featSel = c(featSel1, featSel2, featSel3)
    resultsData[[i]] = list(expVar=expVar, featSel=featSel)
    
  }
  return(resultsData)
}

saveRDS(file="1_results_simulation_data_dim_analysis.rds", resultDataRepeats)







########################
#2, collect resuslts
########################
expVarMat=list()
featSelMat=list()
resultsData=readRDS(file="1_results_simulation_data_dim_analysis.rds")
for (i in 1:length(resultsData)){
  resultRep = resultsData[[i]]
  expVarMat[[i]] = sapply(resultRep, function(x)rowMeans(x[[1]]))
  featSelMat[[i]] = sapply(resultRep, function(x)x[[2]])
}
expVarMat = Reduce("+", expVarMat)/10
featSelMat = Reduce("+", featSelMat)/10





########################
#3, plot results
########################
propSubs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
data <- data.frame(expVarMat)
names(data) <- as.character(propSubs)
barplot(height=as.matrix(data), main="The comparison of prediction performance", 
        ylab="explained variance (R2)", beside=TRUE, col=rainbow(8), xlab="The ratio: subjects number/feature number")
legend("topleft", c("MTLComb-feat","MTLBin", "ridge","lasso", "randomforest", "svm-linear", "svm-radial", "svm-polynomial"), cex=1.0, bty="n", fill=rainbow(8))

data <- data.frame(featSelMat)
names(data) <- as.character(propSubs)
barplot(height=as.matrix(data), main="The comparison of feature selection accuracy", 
        ylab="feature selection accuracy", beside=TRUE, col=rainbow(8), xlab="The ratio: subjects number/feature number")
legend("topleft", c("MTLComb-feat","MTLBin", "ridge + meta","lasso + meta", "randomforest + meta", "svm-linear + meta", "svm-radial + meta", "svm-polynomial + meta"), cex=1.0, bty="n", fill=rainbow(8))












