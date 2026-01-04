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
library(tidyverse)



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
#t=10; p=200; n=100; ctasks=1:5; noise=0.1; props=21
simulateData4MTLComb = function(t=10, p=200, n=100, props=21, ctasks=1:5, noise=0.1, imbRate=0.5){
  
  rtasks=setdiff(1:t, ctasks)
  W <- matrix(data=sign(rnorm(t*p))*rnorm(t*p,mean=1),ncol=t, nrow = p)
  W[props:p,] <- 0
  
  X <- list(); Y <- list(); tX <- list(); tY <- list()
  
  for(i in 1:t){
    X[[i]]  <- matrix(rnorm(n*p), ncol=p, nrow=n)
    tX[[i]] <- matrix(rnorm(n*p), ncol=p, nrow=n)
    
    if (i %in% ctasks){      
      s  <- drop(X[[i]]  %*% W[, i] + noise * rnorm(n))
      ts <- drop(tX[[i]] %*% W[, i] + noise * rnorm(n))
      
      if (imbRate <= 0) {
        Y[[i]]  <- rep(0L, n); tY[[i]] <- rep(0L, n)
      } else{
        if (imbRate >= 0.5) imbRate=0.5
        thr  <- as.numeric(quantile(s,  probs=1-imbRate, type=7))
        tthr <- as.numeric(quantile(ts, probs=1-imbRate, type=7))
        Y[[i]]  <- (sign(s  > thr)-0.5)*2
        tY[[i]] <- (sign(ts > tthr)-0.5)*2
      }
    }else {
      Y[[i]] <- X[[i]] %*% W[,i] + noise * rnorm(n)
      tY[[i]] <- tX[[i]] %*% W[, i] + rnorm(n) * noise
    }
  }
  X=lapply(X, function(x)scale(x))
  tX=lapply(tX, function(x)scale(x))
  Y[rtasks]=lapply(Y[rtasks], function(x)scale(x))
  Y[rtasks]=lapply(Y[rtasks], function(x)scale(x))
  simuData = list(X=X, Y=Y, tX=tX, tY=tY, W=W, ctasks=ctasks, imbRate=imbRate)
  return(simuData)
}

########################
#1, run analysis
########################
opts <- list(init=0, tol=10^-3, maxIter=500, ter=2)
imbRate = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
nt=c(4, 8, 12, 16, 20)

resultDataRepeats=c()
for (rep in 1:10){
  resultData = c()
  for (j in 1:length(nt)){
    print(j)
    rImb = foreach (i = 1:length(imbRate))%dopar%{
      #set.seed(123456)
      rData = simulateData4MTLComb(t=nt[j], p=500, n=100, props=51, ctasks=1:(nt[j]/2), noise=0.5, imbRate=imbRate[i])
      X=rData$X
      Y=rData$Y
      tX=rData$tX
      tY=rData$tY
      ntasks=length(X)
      ctasks=rData$ctasks
      rtasks=setdiff(1:ntasks, ctasks)
      numSigFeat = sum(rData$W[,1]!=0)
      sigFeatIdx = which(rData$W[,1]!=0)
      
      #############
      #MTLCombFeat
      #############
      cvResult = MTLComb_CV(X=X, Y=Y, nfolds=5, lam_ratio=0.01, nlambda=200, ctasks=ctasks, opts=opts, C2=0.1, stratify=TRUE)
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
      bY[rtasks] = lapply(bY[rtasks], function(x)matrix(sign(x-quantile(x,  probs=1-imbRate[i], type=7))))
      cvResult = MTL_L21_CVInSite(X=X, Y=bY, type="classify", nfolds=5, lam_ratio=0.01, nlambda=200, opts=opts, stratify = TRUE)
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
      
      expVar = rbind(expVar1, expVar2)
      featSel = c(featSel1, featSel2)
      return(list(expVar=expVar, featSel=featSel))
    }
    names(rImb)=imbRate
    resultData[[j]] = rImb
  }
  names(resultData) = nt
  resultDataRepeats[[rep]] = resultData
}

saveRDS(file="2_results_simulation_data_imbalance_analysis.rds", resultDataRepeats)




########################
#2, collect resuslts
########################
resultDataRepeats=readRDS(file="2_2_results_simulation_data_imbalance_analysis.rds")

expVarMatRep = vector("list", length(resultDataRepeats))
featSelMatRep = vector("list", length(resultDataRepeats))
for (i in 1:length(resultDataRepeats)){
  resultRep = resultDataRepeats[[i]]
  expVarMat <- vector("list", length(nt))
  featSelMat <- vector("list", length(nt))
  for (j in 1:length(nt)){
    rImb = resultRep[[j]]
    expVarMat[[j]] = sapply(rImb, function(x)rowMeans(x[[1]]))
    featSelMat[[j]] = sapply(rImb, function(x)x[[2]])
  }
  expVarMatRep[[i]] = expVarMat
  featSelMatRep[[i]] = featSelMat
}

expVarMat=vector("list", length(nt))
featSelMat=vector("list", length(nt))
for (j in 1:length(nt)){
  expVarMat[[j]] = Reduce("+", sapply(expVarMatRep, function(x)x[j]))/10
  featSelMat[[j]] = Reduce("+", sapply(featSelMatRep, function(x)x[j]))/10
}
names(featSelMat)=nt
names(expVarMat)=nt







########################
#3, plot results
########################
imb <- c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)

df <- map_dfr(names(featSelMat), function(nt){
  mat <- featSelMat[[nt]]
  tibble(
    nt = as.factor(nt),
    imbalance = imb,
    MTLComb = mat[1,],
    MTLBin  = mat[2,]
  ) %>%
    pivot_longer(
      cols = c(MTLComb, MTLBin),
      names_to = "Method",
      values_to = "Accuracy"
    )
})

p1 <-ggplot(df, aes(x = imbalance, y = Accuracy,
               color = Method, linetype = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  ## Facet titles: rename to "Task number: X"
  facet_wrap(
    ~ nt, nrow = 1,
    labeller = labeller(nt = function(x) paste("Task number:", x))
  ) +
  
  ## Manually set colors and line types
  scale_color_manual(
    values = c(
      "MTLComb" = "blue",
      "MTLBin"  = "red"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "MTLComb" = "solid",
      "MTLBin"  = "dashed"
    )
  ) +
  
  scale_y_continuous(limits = c(0, 1)) +
  
  labs(
    x = "Imbalance Ratio",
    y = "Feature Selection Accuracy",
    color = "Method",
    linetype = "Method"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
ggsave("2_feature_selection_accuracy.svg", p1,
       width = 10, height = 3, device = "svg")




imb <- c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)

df_exp <- map_dfr(names(expVarMat), function(nt){
  mat <- expVarMat[[nt]]
  tibble(
    nt = as.factor(nt),
    imbalance = imb,
    MTLComb = as.numeric(mat["expVar1", ]),  # expVar1: first row (our method)
    MTLBin  = as.numeric(mat["expVar2", ])   # expVar2: second row (baseline method)
  ) %>%
    pivot_longer(
      cols = c(MTLComb, MTLBin),
      names_to = "Method",
      values_to = "ExplainedVar"
    )
})

p2 <-ggplot(df_exp, aes(x = imbalance, y = ExplainedVar,
                   color = Method, linetype = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  facet_wrap(
    ~ nt, nrow = 1,
    labeller = labeller(nt = function(x) paste("Task number:", x))
  ) +
  
  scale_color_manual(
    values = c("MTLComb" = "blue", "MTLBin" = "red")
  ) +
  scale_linetype_manual(
    values = c("MTLComb" = "solid", "MTLBin" = "dashed")
  ) +
  
  # Explained variance can theoretically be negative.
  # Values of expVar2 at small imbalance ratios are close to zero (~1e-15),
  # so we do not enforce [0, 1] limits to avoid truncation.
  # Uncomment the following line if you want to restrict the range to [0, 1].
  # scale_y_continuous(limits = c(0, 1)) +
  
  labs(
    x = "Imbalance Ratio",
    y = "Explained Variance (Prediction Accuracy)",
    color = "Method",
    linetype = "Method"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
ggsave("2_explained_variance.svg", p2,
       width = 10, height = 3, device = "svg")

