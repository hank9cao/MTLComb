rm(list=ls())
source('MTComb_L21.R')

#load data
load('simulated_data.rda')

#z-score normalization
Y[1:10] <- lapply(Y[1:10],
                  function(x) scale(x, center = TRUE, scale = TRUE))
tY[1:10] <- lapply(tY[1:10],
                  function(x) scale(x, center = TRUE, scale = TRUE))
X <- lapply(X, function(x)
    apply(x, 2, function(xx) scale(xx, center = TRUE, scale = TRUE)))
tX <- lapply(tX, function(x)
    apply(x, 2, function(xx) scale(xx, center = TRUE, scale = TRUE)))

#parameters
str(Y)
opts <- list(init=0, tol=10^-6, maxIter=10000, ter=2)
ctasks=ctasks
rtasks = setdiff(1:length(X), ctasks)


####################
#test solver
####################
s1 = solve_L21(X, Y, lam=0.5, C=0, C2=0, ctasks=ctasks, opts=opts)
sum(sqrt(rowSums(s1$W[,ctasks]^2))!=0) #20
sum(sqrt(rowSums(s1$W[,rtasks]^2))!=0) #20


s2 = solve_L21 (X, Y, lam=0.5, C=0, C2=10, ctasks=ctasks, opts=opts)
sum(sqrt(rowSums(s2$W[,ctasks]^2))!=0) #37
sum(sqrt(rowSums(s2$W[,rtasks]^2))!=0) #37

s3 = solve_L21 (X, Y, lam=1, C=0, C2=0, ctasks=ctasks, opts=opts)
sum(sqrt(rowSums(s3$W[,ctasks]^2))!=0) #6
sum(sqrt(rowSums(s3$W[,rtasks]^2))!=0) #6


s4 = solve_L21 (X, Y, lam=0.5, C=100, C2=0, ctasks=ctasks, opts=opts)
sum(cor(s1$W)<0)
#180
sum(cor(s4$W)<0)
#110

plot(s1$Obj, xlab="interations", ylab="objective value")
plot(s2$Obj, xlab="interations", ylab="objective value")
plot(s3$Obj, xlab="interations", ylab="objective value")
plot(s4$Obj, xlab="interations", ylab="objective value")




####################
#test training
####################
#use case 1
fit1=MTLComb_Train(X=X, Y=Y, nlambda=50, lam_ratio=0.1, C=1, C2=0.01, ctasks=ctasks, opts=opts)
#plot regularization tree
matplot(t(sapply(fit1$ws, function(x)sqrt(rowSums(x^2)))), type="l", xlab="lambda sequence", ylab="importance of features")

#show the number of selected  feature for each lambda
plot(sapply(fit1$ws, function(x)sum(sqrt(rowSums(x[,ctasks]^2))!=0)), type = "l", xlab = "lambdas index", 
     ylab = "the number of selected features")
lines(sapply(fit1$ws, function(x)sum(sqrt(rowSums(x[,rtasks]^2))!=0)), col="red")
#show the penalty value for each lambda
plot(sapply(fit1$ws, function(x)sum(sqrt(rowSums(x[,ctasks]^2)))), type = "l", xlab = "lambdas index", 
     ylab = "the penalty value")
lines(sapply(fit1$ws, function(x)sum(sqrt(rowSums(x[,rtasks]^2)))), col="red")

#show the solution of lam_max
sum(fit1$ws[[1]]!=0) #0
rowSums(fit1$ws[[2]]!=0)


#use case 2
fit2=MTLComb_Train(X=X, Y=Y, lambda=0.1, C=1, C2=0.01,ctasks=ctasks, opts=opts)
plot(fit2$Obj, xlab="interations", ylab="objective value")
image(fit2$ws[[1]], xlab="features", ylab="tasks")


#use case 3
fit3=MTLComb_Train(X=X, Y=Y, lambda=c(1, 0.1, 0.01, 0.001), C=1, ctasks=ctasks, opts=opts)
matplot(t(sapply(fit3$ws, function(x)sqrt(rowSums(x^2)))), type="l", xlab="lambda sequence", ylab="importance of features")

#intercept model
fit4=MTLComb_Train(X=X, Y=Y, nlambda=50, lam_ratio=0.001, ctasks=ctasks, C=1, C2=0.1, opts=opts, intercept = T)
dim(fit4$ws[[40]])
#intercepts
fit4$ws[[40]][1,]



#####################
#test cross-validation
#####################
cvResult=MTLComb_CV(X=X, Y=Y, nfolds=10, lam_ratio=0.01, C=1, C2=0.1, nlambda=20, ctasks=ctasks, opts=opts)
plot_MTLComb_CV(cvResult)
cvResult$lambda.regress.min
#0.3095806
cvResult$lambda.classify.min
#0.5026898
cvResult$lambda.weighted.min
#0.4061352

fit5=MTLComb_Train(X=X, Y=Y, nlambda=50, lambda=cvResult$lam_seq[5], ctasks=ctasks, opts=opts)
fit6=MTLComb_Train(X=X, Y=Y, nlambda=50, lambda=cvResult$lambda.weighted.min, ctasks=ctasks, opts=opts)
par(mfrow=c(2,1))
image(fit5$ws[[1]], xlab="features", ylab="tasks", main="lambda selected mannuelly")
image(fit6$ws[[1]], xlab="features", ylab="tasks", main="lambda selected automatically")
par(mfrow=c(1,1))



#####################
#test prediction
#####################
str(tY)
predScores=predict_MTLComb(fit5, newx=tX[[5]], type="regress")
plot(rowMeans(predScores[[1]]), tY[[5]], xlab="averaged prediction scores", ylab="true outcome")
predScores=predict_MTLComb(fit5, newx=tX[[15]], type="classify")
plot(as.factor(tY[[15]]), rowMeans(predScores[[1]]), xlab="classes", ylab="averaged prediction probability")


