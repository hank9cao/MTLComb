LS_MTL_L21 <- function (X, Y, lam, C, opts){
  #min_W sum_k ||Y_k - X_k*W_k||_2^2 +lam||W||_21 + C||W||^2_2 
  #with starting points: 0 or W_0; Data: X and Y; Hyper-parameters: lam, c
  
  proximal_l21 <- function (W, lambda ){
    #prox_f(x)=(1-lam/(max{||x||, lam}))x
    Fnorm <- sqrt(rowSums(W^2))
    zeros <- which(Fnorm==0)              
    prox_f_coef <- 1 - lambda/Fnorm
    prox_f_coef[zeros] <- 0
    prox_f_coef <- prox_f_coef *( prox_f_coef>0 )
    Wp = matrix(rep(prox_f_coef, ncol(W)), nrow=length(prox_f_coef))*W
    return(Wp)
  }
  nonsmooth_eval <- function (W){
    return(sum(sqrt(rowSums(W^2)))*lam)
  }  
  LS_iter_update_tasks <- function(W){
    grad_Ws=sapply(1:nTasks, function(k){(XtX[[k]] %*% W[,k]  - XtY[[k]])/nSubs[k]}) + 2* C * W
    funcVals=LS_funcVal_eval_tasks(W)
    return(list(grad_Ws, funcVals))
  }
  LS_funcVal_eval_tasks <- function ( W){
    funcVal <- sum(sapply(1:nTasks, function(k){mean(0.5 * (Y[[k]] - X[[k]] %*% W[,k])^2)})) + C * norm(W, 'f')^2
    return(funcVal)
  }
  
  #################################  
  # Main algorithm
  #################################  
  Obj <- vector(); 
  nFeats=ncol(X[[1]])
  nTasks=length(X)
  nSubs=sapply(X, nrow)
  log.niterCall=0; log.nfuncCall=0
  XtX=lapply(1:nTasks, function(k)t(X[[k]])%*%X[[k]])
  XtY=lapply(1:nTasks, function(k)t(X[[k]])%*%Y[[k]])

  #initialize a starting point
  if(opts$init==0){
    w0=matrix(0,nrow=nFeats, ncol=length(X))
  }else if(opts$init==1){
    w0 <- opts$w0
  }    
  
  bFlag <- 0; 
  wz <- w0;
  wz_old <- w0;

  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;

  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    ws <- (1 + alpha) * wz - alpha * wz_old;

    # compute function value and gradients of the search point
    iter_update=LS_iter_update_tasks(ws)  
    Gws <- iter_update[[1]]
    Fs <- iter_update[[2]]
    log.niterCall=log.niterCall+1
  
    # the Armijo Goldstein line search scheme
    while (TRUE){
      wzp <- proximal_l21(ws - Gws/gamma, lam / gamma);
      Fzp=LS_funcVal_eval_tasks(wzp)
      log.nfuncCall=log.nfuncCall+1
    

      delta_wzp <- wzp - ws;
      r_sum <- norm(delta_wzp, 'f')^2

      #second order approximation
      Fzp_gamma = Fs + sum(delta_wzp * Gws) + gamma * r_sum/2
    
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
    }
  
    wz_old = wz; wz = wzp; Obj = c(Obj, Fzp + nonsmooth_eval(wz));
  
    #test stop condition.
    if (bFlag) break;
      if (iter>=2){
      if (opts$ter==1 & abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol){
        break
      } else if(opts$ter==2 & abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol*Obj[length(Obj)-1]){
        break
      } else if(opts$ter==3 & iter==opts$maxIter){
        break
      }
    }
  
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
  }
  return(list(W=wzp, Obj=Obj, Logs=c(log.niterCall, log.nfuncCall), gamma=gamma))
}





LR_MTL_L21 <- function (X, Y, lam, C, opts){
  #min_W sum_k{log(1+exp(-Y_k*X_k*W_k))} +lam||W||_21 + c||W||^2_2 
  #with starting points: 0 or W_0; Data: X and Y; Hyper-parameters: lam, C
  
  proximal_l21 <- function (W, lambda ){
    #prox_f(x)=(1-lam/(max{||x||, lam}))x
    Fnorm <- sqrt(rowSums(W^2))
    zeros <- which(Fnorm==0)              
    prox_f_coef <- 1 - lambda/Fnorm
    prox_f_coef[zeros] <- 0
    prox_f_coef <- prox_f_coef *( prox_f_coef>0 )
    Wp = matrix(rep(prox_f_coef, ncol(W)), nrow=length(prox_f_coef))*W
    return(Wp)
  }
  nonsmooth_eval <- function (W){
    return(sum(sqrt(rowSums(W^2)))*lam)
  }  
  LR_iter_update_tasks <- function(W){
    rData=lapply(1:nTasks, function(k){
      x=X[[k]]; y=Y[[k]];w=W[,k]
      weight <- 1/length(y)
      l <- -y*(x %*% w) 
      lp <- l
      lp[lp<0] <- 0
      funcVal <- sum(weight * ( log( exp(-lp) +  exp(l-lp) ) + lp ))
      b <- (-weight*y)*(1 - 1/ (1+exp(l)))
      grad_w <- (t(x) %*% b)
      return(list(grad_w=grad_w, funcVal=funcVal))
    })
    grad_Ws=sapply(rData, function(x)x$grad_w) + 2* C * W
    funcVals=sum(sapply(rData, function(x)x$funcVal)) + C * norm(W, 'f')^2
    return(list(grad_Ws=grad_Ws, funcVals=funcVals))
  }
  
  LR_funcVal_eval_tasks <- function (W){
    funcVal <- sum(sapply(1:nTasks, function(k){
      x=X[[k]]; y=Y[[k]];w=W[,k]
      weight <- 1/length(y)
      l <- - y*(x %*% w)
      lp <- l
      lp[lp<0] <- 0
      return(sum(weight * ( log( exp(-lp) +  exp(l-lp) ) + lp )))
    })) + C * norm(W, 'f')^2
    return(funcVal)
  }

  #################################  
  # Main algorithm
  #################################  
  Obj <- vector(); 
  nFeats=ncol(X[[1]])
  nTasks=length(X)
  log.niterCall=0; log.nfuncCall=0
  
  #initialize a starting point
  if(opts$init==0){
    w0=matrix(0,nrow=nFeats, ncol=length(X))
  }else if(opts$init==1){
    w0 <- opts$w0
  }    
  
  bFlag <- 0; 
  wz <- w0;
  wz_old <- w0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    ws <- (1 + alpha) * wz - alpha * wz_old;
    
    # compute function value and gradients of the search point
    iter_update=LR_iter_update_tasks(ws)  
    Gws <- iter_update$grad_Ws
    Fs <- iter_update$funcVals
    log.niterCall=log.niterCall+1
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      wzp <- proximal_l21(ws - Gws/gamma, lam / gamma);
      Fzp=LR_funcVal_eval_tasks(wzp)
      log.nfuncCall=log.nfuncCall+1
      
      delta_wzp <- wzp - ws;
      r_sum <- norm(delta_wzp, 'f')^2
      
      #second order approximation
      Fzp_gamma = Fs + sum(delta_wzp * Gws) + gamma * r_sum/2
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
    }
    
    wz_old = wz; wz = wzp; Obj = c(Obj, Fzp + nonsmooth_eval(wz));
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (opts$ter==1 & abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol){
        break
      } else if(opts$ter==2 & abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol*Obj[length(Obj)-1]){
        break
      } else if(opts$ter==3 & iter==opts$maxIter){
        break
      } 
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
  }
  return(list(W=wzp, Obj=Obj, Logs=c(log.niterCall, log.nfuncCall), gamma=gamma))
}






MTL_L21_Train = function(X=NULL, Y=NULL, type="regress", nlambda=10, lam_ratio=0.01, lambda=NULL, C=0, 
                         opts=list(init=0, maxIter=20, tol=0.01, ter=2)){
  #initialize final result
  fit=list();fit$ws=list();fit$Logs=vector();fit$Obj=vector();fit$gamma=vector();fit$type=type
  nTasks=length(X)
  xys=sapply(1:nTasks, function(k)t(X[[k]])%*%Y[[k]]/nrow(X[[k]]))
  xy_norm=sqrt(rowSums(xys^2))
  
  #########################
  #for regression
  #########################
  if (type=="regress"){
    #########################
    #Train the solution path of lasso
    #########################
    #estimate lambda sequence
    if(length(lambda)>1){
      lam_seq=lambda
    } else if(length(lambda)==1){
      lam_max=max(xy_norm)
      lam_min=lambda
      lam_seq=exp(seq(log(lam_max),log(lam_min),length.out = nlambda))
    } else if(is.null(lambda)){
      lam_max=max(xy_norm)
      lam_min=lam_ratio*lam_max
      lam_seq=exp(seq(log(lam_max),log(lam_min),length.out = nlambda))
    }

    #warm-start training procedure
    optsTrain=opts
    for(i in 1:length(lam_seq)){
      m=LS_MTL_L21(X=X, Y=Y, lam=lam_seq[i], C=C, opts=optsTrain)
      optsTrain$w0=m$W; optsTrain$init=1
      fit$ws[[i]]=m$W
      fit$Obj=c(fit$Obj, m$Obj)
      fit$Logs=rbind(fit$Logs, m$Logs)
      fit$gamma=c(fit$gamma, m$gamma)
    }
    fit$lam_seq=lam_seq
    names(fit$ws)=paste0("Lam=", round(lam_seq, 2))
  }
  #########################
  #for classification
  #########################
  else if(type=="classify"){
    #########################
    #Train the solution path of lasso 
    #########################
    #estimate lambda sequence
    if(length(lambda)>1){
      lam_seq=lambda
    } else if(length(lambda)==1){
      lam_max=max(xy_norm/2)
      lam_min=lambda
      lam_seq=exp(seq(log(lam_max),log(lam_min),length.out = nlambda))
    } else if(is.null(lambda)){
      lam_max=max(xy_norm/2)
      lam_min=lam_ratio*lam_max
      lam_seq=exp(seq(log(lam_max),log(lam_min),length.out = nlambda))
    }
    
    #warm-start training procedure
    optsTrain=opts
    for(i in 1:length(lam_seq)){
      m=LR_MTL_L21(X=X, Y=Y, lam=lam_seq[i], C=C, opts=optsTrain)
      optsTrain$w0=m$W; optsTrain$init=1
      fit$ws[[i]]=m$W
      fit$Obj=c(fit$Obj, m$Obj)
      fit$Logs=rbind(fit$Logs, m$Logs)
      fit$gamma=c(fit$gamma, m$gamma)
    }
    
    fit$lam_seq=lam_seq
    names(fit$ws)=paste0("Lam=", round(lam_seq, 2))
  }
  
  fit$ws=lapply(fit$ws, function(x){rownames(x)=colnames(X[[1]]); return(x)})
  if(length(lambda)==1){fit$ws=fit$ws[length(fit$ws)]}

  return(fit)
}




MTL_L21_CVCroSite = function(X=NULL, Y=NULL, type="regress", lam_ratio=0.01, nlambda=10, lambda=NULL,
                             opts=list(init=0, maxIter=50, tol=0.001, ter=2), C=0){
  nTasks=length(X)
  cvResult=list(); cvResult$type=type; cvResult$C=C; 

  if (type=="regress"){
    mse_fold=vector()
    lam_seq=vector()
    for (k in 1:nTasks){

      Xtest=X[[k]]; Ytest=Y[[k]]
      Xtrain=X[-k]; Ytrain=Y[-k]
      
      fit=MTL_L21_Train(X=Xtrain, Y=Ytrain, nlambda=nlambda, lam_ratio=lam_ratio, type="regress", opts=opts, C=C, lambda=lambda)
      yhat=sapply(fit$ws, function(x) rowMeans(Xtest%*%x))
      mse_fold=rbind(mse_fold, colMeans(apply(yhat, 2, function(x)(x-Ytest)^2)))
      lam_seq=rbind(lam_seq, fit$lam_seq)
    }
    lambda.min=mean(sapply(1:nTasks, function(k)lam_seq[k,order(mse_fold[k,])[1]]))
    colnames(mse_fold)=paste0("Lam",1:ncol(mse_fold))
    colnames(lam_seq)=paste0("Lam",1:ncol(mse_fold))
    cvResult$lam_seq=lam_seq; cvResult$mse_fold=mse_fold; cvResult$lambda.min=lambda.min
  } else if(type=="classify"){
    mcr_fold=vector()
    lam_seq=vector()
    for (k in 1:nTasks){
      
      Xtest=X[[k]]; Ytest=Y[[k]]
      Xtrain=X[-k]; Ytrain=Y[-k]
      
      fit=MTL_L21_Train(X=Xtrain, Y=Ytrain, nlambda=nlambda, lam_ratio=lam_ratio, type="classify", opts=opts, C=C, lambda=lambda)
      yhat=sapply(fit$ws, function(x) rowMeans(1/(1+exp(-Xtest%*%x))))
      yhat=sign(yhat-0.5); yhat[yhat==0]=1
      mcr_fold=rbind(mcr_fold, apply(yhat, 2, function(x)mean(x!=Ytest)))
      lam_seq=rbind(lam_seq, fit$lam_seq)
    }
    lambda.min=mean(sapply(1:nTasks, function(k)lam_seq[k,order(mcr_fold[k,])[1]]))
    colnames(mcr_fold)=paste0("Lam",1:ncol(mcr_fold))
    colnames(lam_seq)=paste0("Lam",1:ncol(mcr_fold))
    cvResult$lam_seq=lam_seq; cvResult$mcr_fold=mcr_fold; cvResult$lambda.min=lambda.min
  }
  return(cvResult)
}


getCVPartition <- function(Y, cv_fold, stratify){
  task_num = length(Y);
  randIdx <- lapply(Y, function(x) sample(1:length(x), length(x), replace = FALSE))   
  
  cvPar = {};
  for (cv_idx in 1: cv_fold){
    # build cross validation data splittings for each task.
    cvTrain = {};
    cvTest = {};
    
    #stratified cross validation
    for (t in 1: task_num){
      task_sample_size <- length(Y[[t]]);
      if (stratify){
        ct <- which(Y[[t]][randIdx[[t]]]==-1);
        cs <- which(Y[[t]][randIdx[[t]]]==1);
        ct_idx <- seq(cv_idx, length(ct), cv_fold);
        cs_idx <- seq(cv_idx, length(cs), cv_fold);
        te_idx <- c(ct[ct_idx], cs[cs_idx]);
        tr_idx <- seq(1,task_sample_size)[!is.element(1:task_sample_size, te_idx)];
      } else {
        te_idx <- seq(cv_idx, task_sample_size, by=cv_fold)
        tr_idx <- seq(1,task_sample_size)[!is.element(1:task_sample_size, te_idx)]
      }
      cvTrain[[t]] = randIdx[[t]][tr_idx]
      cvTest[t] = list(randIdx[[t]][te_idx])
    }
    cvPar[[cv_idx]]=list(cvTrain=cvTrain, cvTest=cvTest);
  }
  return(cvPar)
}

MTL_L21_CVInSite = function(X=NULL, Y=NULL, type="regress", nfolds=10, stratify=F, lam_ratio=0.01, nlambda=10, lambda=NULL,
                             opts=list(init=0, maxIter=50, tol=0.001, ter=2), C=0){
  nTasks=length(X)
  cvResult=list(); cvResult$type=type; cvResult$C=C; 
  cvPar <- getCVPartition(Y, nfolds, stratify)
  
  if (type=="regress"){
    mse_fold=vector()
    lam_seq=vector()
    for (i in 1:length(cvPar)){
      Xtrain=lapply(c(1:nTasks), function(k) X[[k]][cvPar[[i]][[1]][[k]], ])
      Ytrain=lapply(c(1:nTasks), function(k) Y[[k]][cvPar[[i]][[1]][[k]]])
      Xtest=lapply(c(1:nTasks), function(k) X[[k]][cvPar[[i]][[2]][[k]], ])
      Ytest=lapply(c(1:nTasks), function(k) Y[[k]][cvPar[[i]][[2]][[k]]])
      
      fit=MTL_L21_Train(X=Xtrain, Y=Ytrain, nlambda=nlambda, lam_ratio=lam_ratio, type="regress", opts=opts, C=C, lambda=lambda)
      yhat=lapply(fit$ws, function(w)lapply(1:nTasks, function(x) Xtest[[x]]%*%w[,x]))
      mse_fold=rbind(mse_fold, sapply(yhat, function(x)mean(sapply(1:nTasks, function(xx)mean((x[[xx]]-Ytest[[xx]])^2)))))
      lam_seq=rbind(lam_seq, fit$lam_seq)
    }
    lambda.min=colMeans(lam_seq)[order(colMeans(mse_fold))[1]]
    colnames(mse_fold)=paste0("Lam",1:ncol(mse_fold))
    colnames(lam_seq)=paste0("Lam",1:ncol(mse_fold))
    cvResult$lam_seq=lam_seq; cvResult$mse_fold=mse_fold; cvResult$lambda.min=lambda.min
    
  } else if(type=="classify"){
    mcr_fold=vector()
    lam_seq=vector()
    for (i in 1:length(cvPar)){
      
      Xtrain=lapply(c(1:nTasks), function(k) X[[k]][cvPar[[i]][[1]][[k]], ])
      Ytrain=lapply(c(1:nTasks), function(k) Y[[k]][cvPar[[i]][[1]][[k]]])
      Xtest=lapply(c(1:nTasks), function(k) X[[k]][cvPar[[i]][[2]][[k]], ])
      Ytest=lapply(c(1:nTasks), function(k) Y[[k]][cvPar[[i]][[2]][[k]]])
      
      fit=MTL_L21_Train(X=Xtrain, Y=Ytrain, nlambda=nlambda, lam_ratio=lam_ratio, type="classify", opts=opts, C=C, lambda=lambda)
      yhat=lapply(fit$ws, function(w)lapply(1:nTasks, function(x) {lab=sign(Xtest[[x]]%*%w[,x]); lab[lab==0]=1; return(lab)}))
      mcr_fold=rbind(mcr_fold, sapply(yhat, function(x)mean(sapply(1:nTasks, function(xx)mean(x[[xx]]!=Ytest[[xx]])))))
      lam_seq=rbind(lam_seq, fit$lam_seq)
    }
    lambda.min=colMeans(lam_seq)[order(colMeans(mcr_fold))[1]]
    colnames(mcr_fold)=paste0("Lam",1:ncol(mcr_fold))
    colnames(lam_seq)=paste0("Lam",1:ncol(mcr_fold))
    cvResult$lam_seq=lam_seq; cvResult$mcr_fold=mcr_fold; cvResult$lambda.min=lambda.min
  }
  return(cvResult)
}



MTL_L21_predict=function(fit, newx){
  if(fit$type=="regress"){
    yhat=lapply(fit$ws, function(x)newx%*%x)
    
  }else if(fit$type=="classify"){
    yhat=lapply(fit$ws, function(x)1/(1+exp(-newx%*%x)))
    
  }
  return(yhat)
}
