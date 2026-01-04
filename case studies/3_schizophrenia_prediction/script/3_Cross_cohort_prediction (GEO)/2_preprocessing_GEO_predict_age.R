rm(list=ls())
gc()


#0, remove subjects with ages>65 and <18
#1, remove outliers 
#2, combine
#3, sva analysis
#4, get 5 pc
#5, remove covariates
#6, seperate


load("../../preprocessed_data/GEO/unnormalized_X_Y.rda")

#construct the covariate matrix
for(i in 1:3){
  Y[[i]]$pheno <- as.factor(Y[[i]]$pheno)  
  Y[[i]]$pmi <- as.numeric(Y[[i]]$pmi)  
  Y[[i]]$gender <- as.factor(as.character((Y[[i]]$gender)))
  X[[i]] = X[[i]][Y[[i]]$ages<=65, ]
  Y[[i]] = Y[[i]][Y[[i]]$ages<=65, ]
  X[[i]] = X[[i]][Y[[i]]$pheno==0, ]
  Y[[i]] = Y[[i]][Y[[i]]$pheno==0, ]
}

#remove outliers
for (i in 1:3){
  curPCA = prcomp(X[[i]], scale=TRUE)$x[,1:2]
  out1 = which( abs((curPCA[,1] - mean(curPCA[,1]))/sd(curPCA[,1])) > 4 )
  out2 = which( abs((curPCA[,2] - mean(curPCA[,2]))/sd(curPCA[,2])) > 4 )
  outlier = unique(sort(c(out1, out2)))
  X11();plot(curPCA);points(curPCA[outlier, ], col="red"); 
  if(length(outlier)>0){
    X[[i]]  <- X[[i]][-outlier, ]
    Y[[i]]  <- Y[[i]][-outlier, ]
  }
}

#sva analysis
library('sva')
covY <- Y
sites <- as.factor(unlist(sapply(1:3, function(x)rep(x, nrow(covY[[x]])))))
covY <- cbind(do.call(rbind, covY), sites)
X <- do.call(rbind, X)
mod = model.matrix(~ages+gender+ph+pmi+sites, data=covY)
mod0 = model.matrix(~gender+ph+pmi+sites+1, data=covY)
#n.sv = num.sv(t(X), mod) 
n.sv=5
SVAsol = sva(dat=t(X), mod=mod, mod0=mod0, n.sv=n.sv)
svs = SVAsol$sv
covY <- cbind(covY, svs)
rownames(covY) <- covY[,1]
covY <- covY[,-1]
covY <- covY[,-1]
colnames(covY) <- c("age", "sex", "ph", "pmi", "rin", "site", paste0("sv", 1:n.sv))



#calculate 5 pcs
ns.pca = 5
pctest = prcomp(X, scale=TRUE)
pcas <- pctest$x[,1:ns.pca]
covY <- cbind(covY, pcas)


#Remove covariates
Xres=apply(X, 2, function(x)lm(x~sex+ph+pmi+site+sv1+sv2+sv3+sv4+sv5+PC1+PC2+PC3+PC4+PC5, data=covY)$res)
Xres=scale(Xres)



saveRDS(file="../../preprocessed_data/GEO/gene_expression_predict_age.rds", Xres)
saveRDS(file="../../preprocessed_data/GEO/y_matrix_predict_age.rds", covY)