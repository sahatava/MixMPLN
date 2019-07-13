library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(cvTools)
library(huge)
library(edgeR)
 
##########################################################
#' MixMPLN function
#'
#' This function receives the sample-taxa count matrix and cluser the samples based on the MixMPLN method.  
#' @param file name of the csv file as input data
#' @param K number of components
#' @param penalty "no_sparsity", "CV", "StARS" , "fixed" , "iterative"
#' @param init "Kmeans"
#' @param out "precision" , "partial" , "adj"
#' @param threshold a value between 0 and 1 which indicates the threshod to generate the adjacency matrix
#' @return
#' @export
MixMPLN_function <-function(X,K, penalty , init  ){

   if (penalty == "no_sparsity"){
    package = "no_sparsity"  
  }
  if (penalty == "CV"){
    package = "glasso"  
    penalty = 1
  }
  if (penalty == "StARS"){
    package = "huge"  
    penalty = 1
  }
  if (penalty == "fixed"){
    package = "huge"  
    penalty = 2
  }
  if (penalty == "iterative"){
    package = "huge"  
    penalty = 3
  }
  #X = read.csv(file)
  #X = as.matrix(X[, 2:dim(X)[2]])
  X = X+1
  niter=100
  diff_threshold = .01
  if(K<=0){
    stop ("Error: K, number of components, needs to be an integer >= 1")
  }
  N=dim(X)[1]
  M=dim(X)[2]
  
  if((N/K)<=2){
    stop("Error: too many clusters and not enough data points!")
  }
  
  pi_j=rep(1/K,K)

  sigma=array(rep(0,K*M*M),c(K,M,M))
  lambda=array(rep(0,N*K*M),c(N,K,M))
  mu=array(rep(0,K*M),c(K,M))
  
  if(init == "Kmeans"){
    member = kmeans(X,K)$cluster
    while(min(table(member))<2){
      ## random assignment when kmeans has a cluster with less than two member
      member=sample.int(K,size=N,replace=TRUE)
    }
  }# if(init = "Kmeans")
  else{
    member = sample.int(K,size=N,replace=TRUE)
    while(min(table(member))<N/K-2){
      ## random assignment when kmeans has a cluster with less than two member
      member=sample.int(K,size=N,replace=TRUE)
    }
  }#else 

    
  
  for(k in 1:K){
    tsigma=matrix(rep(0,M*M),ncol=M)
    ind=which(member==k)
    for(i in 1:M){
      for(j in i:M){
        if(i==j){
          mi=mean(X[ind,i])
          tsigma[i,i]=log(abs( (var(X[ind,i])-mi)/(mi^2+1e-5) ) + 1)
        }else{
          mi=mean(X[ind,i])
          mj=mean(X[ind,j])
          v=log(abs( var(X[ind,i],X[ind,j])/(mi*mj+1e-5) ) + 1)
          tsigma[i,j]=v
          tsigma[j,i]=v
        }
      }
    }
    sigma[k,,]=make.positive.definite(tsigma,tol=1e-5)

    tmu=log(colMeans(X[ind,]))
    mu[k,]=tmu
    lambda[ind,k,1:M]=log(X[ind,1:M]+1)
  }
  
 
  
  invSigma=array(rep(0,K*M*M),c(K,M,M))
  invSigma_init=array(rep(0,K*M*M),c(K,M,M))
  for(k in 1:K){
    invSigma[k,,]=solve(sigma[k,,])
  }
  
  wij=compute_wij(X=X,invSigma=invSigma,lambda=lambda,mu=mu,K=K,PI=pi_j)
  ll=compute_ll_mix(wij)  
  
  norm_wij=compute_norm_wij(wij)
  pi_j=colMeans(norm_wij)   
  for(k in 1:K){
    C=cov.wt(x=lambda[1:N,k,1:M],wt=norm_wij[1:N,k],method="ML")
    invSigma_init[k,,]=solve(make.positive.definite(C$cov,tol=1e-5))
  }
  curr_lambda=lambda
  curr_mu=mu
  curr_invSigma=invSigma
  curr_PI=pi_j
  curr_norm_wij=norm_wij
  curr_ll=ll
  
  if(package=="no_sparsity"){
    invSigma=compute_invSigma_mix(wij=curr_norm_wij,lambda=curr_lambda,K=K)
  }
  ## using the glasso package
  if(package=="glasso"){
    invSigma=compute_invSigma_mix_glasso(wij=curr_norm_wij,lambda=curr_lambda,K=K,invSigma_init,penalty)
  }
  ## using the huge package
  if(package=="huge"){
    invSigma=compute_invSigma_mix_huge(wij=curr_norm_wij,data=curr_lambda,K=K,invSigma_init,penalty)
  }

   
  lambda=compute_lambda_mix(X=X,invSigma=invSigma,mu=curr_mu,K=K,lambda=curr_lambda)
  mu=compute_mu_mix(wij=curr_norm_wij,lambda=lambda)
  
  wij=compute_wij(X=X,invSigma=invSigma,mu=mu,lambda=lambda,K=K,PI=curr_PI)
  lln=compute_ll_mix(wij)
  norm_wij=compute_norm_wij(wij)
  pi_j=colMeans(norm_wij)

  ## 
  diff_sigma = 0
  for(l in 1:K){
    diff_sigma = diff_sigma + myfr(solve(curr_invSigma[l,,]),solve(invSigma[l,,]))
  }
  diff_sigma= diff_sigma/K
  #print(diff_sigma)
  ##
  i=0 
  while(curr_ll<lln && i<niter && diff_sigma>diff_threshold){
    i=i+1
    print(c("iteration" , i))
    curr_lambda=lambda
    curr_mu=mu
    curr_invSigma=invSigma
    curr_PI=pi_j
    curr_norm_wij=norm_wij
    curr_ll=lln

    if(package=="no_sparsity"){
      invSigma=compute_invSigma_mix(wij=curr_norm_wij,lambda=curr_lambda,K=K)
    }
    ## using the glasso package
    if(package=="glasso"){
      invSigma=compute_invSigma_mix_glasso(wij=curr_norm_wij,lambda=curr_lambda,K=K,invSigma_init,penalty)
    }
    ## using the huge package
    if(package=="huge"){
      invSigma=compute_invSigma_mix_huge(wij=curr_norm_wij,data=curr_lambda,K=K,invSigma_init,penalty)
    }

    lambda=compute_lambda_mix(X=X,invSigma=invSigma,mu=curr_mu,K=K,lambda=curr_lambda)
    mu=compute_mu_mix(wij=curr_norm_wij,lambda=lambda)
    wij=compute_wij(X=X,invSigma=invSigma,mu=mu,lambda=lambda,K=K,PI=curr_PI)
    norm_wij=compute_norm_wij(wij)
    pi_j=colMeans(norm_wij)
    lln=compute_ll_mix(wij)
    diff_sigma = 0
    for(l in 1:K){
      diff_sigma = diff_sigma + myfr(solve(curr_invSigma[l,,]),solve(invSigma[l,,]))
    }
    diff_sigma= diff_sigma/K
    #print(c("diff_sigma", diff_sigma))
  }#while
 
 precision_out = list()
 partial_out = list()
 lambda_out = list()

 
   for(l in 1:K){
    precision_out[[l]] = curr_invSigma[l,,] 
    colnames(precision_out[[l]]) = colnames(X)
    partial_out[[l]] = -cov2cor(curr_invSigma[l,,]) 
    colnames(partial_out[[l]]) = colnames(X)
    lambda_out[[l]] = curr_lambda[,l,]
   }
 


  
 lcluster=apply(curr_norm_wij,1,which.max)


return( list("precision"=precision_out , "partial"=partial_out , "lambda"= lambda_out, "pi"=curr_PI , "cluster"=lcluster , "ll" = curr_ll))
}
##########################################################
# -exp(x) - b*x + c
# Newton-Raphson method from numerical recipes in C after modifications to fit the function we have
#' A function to    
#' @param b
#' @param c
#' @param xacc
#' @return x1
#' @export
mysolve<-function(b,c,xacc=1e-3){
   
  if(b<0){
    stop("Error: b cannot be <0")
  }
  rts=0
  #bounding the interval where the root is
  #print(c("c",c))
  if(c==0){
    x1=-800
    x2=0
  }
  if(c>0){
    if(c>=1){
      x1=-log(1+c)
      x2=log(1+c)
    }else{
      x1=log(c)
      x2=log(1+c)
    }
  }else{
    if(c<0){
    x1=-((abs(c)+2)/b)
    x2=0
    }
  }
  
  fl=-exp(x1) - b*x1 + c
  fh=-exp(x2) - b*x2 + c

  
#  if((fl>0 && fh>0)||(fl<0 && fh<0)){
#    stop (c("Error in root finding"," b= ",b," c= ",c))
#  }
  
  if(fl==0){
    return (x1)
  }
  if(fh==0){
    return (x2)
  }

  if(fl<0){
    xl=x1
    xh=x2
  }else{
    xh=x1
    xl=x2
  }
  rts=0.5*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  f=-exp(rts) - b*rts + c
  df=-exp(rts) - b
  
  MAXIT=100
  j=0
  while(j<=MAXIT){
    j=j+1
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0) || (abs(2*f)>abs(dxold*df))){
      dxold=dx
      dx=0.5*(xh-xl)
      rts=xl+dx
      if(xl==rts){
        return (rts)
      }
    }else{
      dxold=dx
      dx=f/df
      temp=rts
      rts=rts - dx
      if(temp==rts){
        return (rts)
      }
    }
    if(abs(dx) < xacc) {
      return (rts)
    }
    f=-exp(rts) - b*rts + c
    df=-exp(rts) - b
    if(f<0){
      xl=rts
    }else{
      xh=rts
    }
  }
  return (x1)
}
#################################################### 
####################################################
#' A function to    
#'    
#'   
#' @param wij
#' @return val
#' @export 
compute_ll_mix<-function(wij){
  N=dim(wij)[1]
  K=dim(wij)[2]
  val=0
  for(i in 1:N){
    x=wij[i,]
    y=max(x)
    for(k in 1:K){
      val = val + y + log(sum(exp(x-y))) 
    }
  }
#  print(c("ll val",val))
  return (val)
}
##################################################
#' A function to  
#' @param X  
#' @param invSigma 
#' @param mu
#' @param lambda
#' @param K
#' @param PI
#' @return wij
#' @export
compute_wij<-function(X,invSigma,mu,lambda,K,PI){
  N=dim(X)[1]
  M=dim(X)[2]
  ldetK=rep(0,K)
  for(k in 1:K){
    ldetK[k]=determinant(invSigma[k,,],logarithm = TRUE)$modulus[1]
  }
  wij=matrix(rep(0,N*K),ncol=K)
  for(i in 1:N){
    for(k in 1:K){
      ld=compute_ldensity(X=X[i,],invSigma=invSigma[k,,],lambda=lambda[i,k,],ldt=ldetK[k],mu=mu[k,])
      if(is.nan(ld)){
        wij[i,k]=0
      }else{
        wij[i,k]= log(PI[k]) + ld
      }
      if(is.nan(ld)){
        print (c("MESS NAN",i," ",k))
      }
    }
  }
  
  
  return(wij)
}
#################################################
#' A function to 
#' @param X 
#' @param invSigma 
#' @param lambda
#' @param ldt
#' @param mu
#' @return val
#' @export
compute_ldensity<-function(X,invSigma,lambda,ldt,mu){
  el = sum(exp(lambda))
  lx = sum(lambda*X)
  mv = (t(lambda-mu)%*%invSigma%*%(lambda-mu))/2
  lf = sum(lfactorial(X))
  val = (ldt/2) -el -lf + lx -mv 
 
  return(val)
}
################################################
#' A function to    
#'    
#'    
#' @param wij
#' @return norm_wij
#' @export
compute_norm_wij<-function(wij){
  N=dim(wij)[1]
  K=dim(wij)[2]
  norm_wij=matrix(rep(0,N*K),ncol=K)
  for(i in 1:N){
    x=wij[i,]
    for(k in 1:K){
      if(is.infinite(x[k])){
        norm_wij[i,k]=0
      }else{
        y=(sum(exp(x-x[k])))
        norm_wij[i,k]=1/(sum(exp(x-x[k])))
      }
    }
  }
  
  for(k in 1:K){
   if(sum(norm_wij[,k]<=1e-100)==N){
      norm_wij[,k] = norm_wij[,k] + 1e-100
    }
  }

 
  norm_wij=t(scale(t(norm_wij),center=FALSE,scale=rowSums(norm_wij)))
  return(norm_wij)
}
#############################################
#' A function to    
#'    
#' @param w_ij
#' @param lambda
#' @return mu
#' @export
compute_mu_mix<-function(wij,lambda){
  N=dim(lambda)[1]
  K=dim(lambda)[2]
  M=dim(lambda)[3]
  mu=array(rep(0,K*M),c(K,M))
  for(k in 1:K){
    for(m in 1:M){
      mu[k,m]=sum(wij[1:N,k]*lambda[1:N,k,m])/sum(wij[1:N,k])
    }
  }
  return(mu)
}
############################################
#' A function to  
#' @param X
#' @param invSigma
#' @param K
#' @param mu
#' @param lambda
#' @return lambda
#' @export
compute_lambda_mix<-function(X,invSigma,K,mu,lambda){
  N=dim(X)[1]
  M=dim(X)[2]

  for(k in 1:K){
    lambdaij=compute_lambdaij_mix(X=X,invSigma=invSigma[k,,],mu=mu[k,],
                                        lambda=lambda[1:N,k,],k=k)
    lambda[1:N,k,]=lambdaij
  }
  return(lambda)
}
#############################################
#' A function to 
#' @param X
#' @param invSigma 
#' @param mu
#' @param lambda 
#' @param K
#' @return lambdanew
#' @export
compute_lambdaij_mix<-function(X,invSigma,mu,lambda,k){
  N=dim(X)[1]
  M=dim(X)[2]
  lambdanew=matrix(rep(0,N*M),ncol=M)
    
  for(i in 1:N){
    S=(lambda[i,]-mu)%*%invSigma
    for(j in 1:M){
      v=lambda[i,j]*invSigma[j,j]
      #print(c("invSigma[j,j]",invSigma[j,j]))
      #print(c("X[i,j]",X[i,j]))
      #print(c("S[j]",S[j]))
      #print(c("v",v))
      lambdanew[i,j]=mysolve(b=invSigma[j,j],c=(X[i,j]-S[j]+v))
 
    }
  }
  return(lambdanew)
}
############################################
#' A function to    
#' @param wij
#' @param lambda
#' @param K
#' @return invSigma
#' @export
compute_invSigma_mix<-function(wij,lambda,K){
  N=dim(lambda)[1]
  M=dim(lambda)[3]
  invSigma=array(rep(0,K*M*M),c(K,M,M))
  for(k in 1:K){
    C=cov.wt(x=lambda[1:N,k,1:M],wt=wij[1:N,k],method="ML")
 
    invSigma[k,,]=solve(make.positive.definite(C$cov,tol=1e-5))
    for (i in 1:M){
            if(is.na(invSigma[k,i,i]) || is.infinite(invSigma[k,i,i]) ){ invSigma[k,i,i]= 0}
    }
    invSigma[k,,] = make.positive.definite(invSigma[k,,] ,tol=1e-5)
 
  }
  return(invSigma)
} 
############################################# 
#' A function to  
#' @param wij
#' @param lambda
#' @param K
#' @param invSigma_init
#' @param penalty
#' @return invSigma
#' @export
compute_invSigma_mix_glasso<-function(wij,lambda,K,invSigma_init,penalty){
  N=dim(lambda)[1]
  M=dim(lambda)[3]
  invSigma=array(rep(0,K*M*M),c(K,M,M))
  lcluster=apply(wij,1,which.max)

  for(k in 1:K){
    id=which(lcluster==k)
    n = length(id)
    C=cov.wt(x=lambda[1:N,k,1:M],wt=wij[1:N,k],method="ML")
    sigma2 = C$cov
    for (i in 1:M){
            if(is.na(sigma2[i,i]) || is.infinite(sigma2[i,i]) ){ sigma2[i,i]= 0}
          }
    Sigma2 = make.positive.definite(sigma2,tol=1e-5)

    ## find the best rho by cross validation
    if(penalty==1){
      rhoMax = glasso_cv_mixture(wij[1:N,k],lambda[1:N,k,1:M],5, rholist=NULL, verbose=T)
      #print(c("rhomax------------",rhoMax ))
    }
    ## find the rho by 2(M)^2/(N||invSigma||)
    if(penalty==2){
      rhoMax = 2*(M^2)/(N*sum(abs(invSigma_init[k,,])))
      #print(c("rhomax------------",rhoMax ))
    }
    ## different rho in each iteration
    if(penalty==3){
      if(n>2){
        rhoMax = 2*(M^2)/(n*sum(abs(solve(make.positive.definite(cov(lambda[id,k,1:M]),tol=1e-5)))))  
      }else{
        rhoMax = 0
      }
      #print(c("rhomax------------",rhoMax ))
    }
    
    GlassoRes = huge( Sigma2 , lambda = rhoMax , method = "glasso", cov.output = TRUE, verbose = TRUE)
    #GlassoRes = glasso( Sigma2 ,rhoMax, start="warm", w.init=Sigma2 , wi.init=solve(Sigma2))
    invSigma[k,,] = as.matrix(GlassoRes$icov[[1]]) 
    for (i in 1:M){
            if(is.na(invSigma[k,i,i]) || is.infinite(invSigma[k,i,i]) ){ invSigma[k,i,i]= 0}
    }
    invSigma[k,,] = make.positive.definite(invSigma[k,,] ,tol=1e-5)
    
  }
  return(invSigma)
}
 
 
################################################################### 
## to find the best penalty in glasso for the mixture model 
#' A function to    
#' @param W
#' @param ts
#' @param k
#' @param rholist 
#' @param verbose 
#' @return make.positive.definite(Sigma_emp ,tol=1e-5)
#' @export
glasso_cv_mixture <- function(W, ts, k , rholist=NULL, verbose=T) {

#if(FALSE){
######### huge
  #low_range = 0.001
  #mu = compute_mu_mixture(W,ts)
  #S = compute_sigma_mixture_emp(W, ts, mu)
  #est = huge(S, nlambda = 10, method = "glasso",lambda.min.ratio=0.005, cov.output = TRUE, verbose = TRUE)
  #rholist = est$lambda
  #print("huge")
  #print(rholist)
#}# if(FALSE)
######### glass
   low_range = 0.001 
 # if (is.null(rholist)) {
    ## We will determine the maximum rho as one where the inverse covariance
    ## has at least 5% non-zero off-diagonal elements (arbitrary)
    mu = compute_mu_mixture(W,ts)
    S = compute_sigma_mixture_emp(W, ts, mu)
    rholist = seq(low_range, max(abs(S))+low_range , length=11)[-1]
    M = dim(ts)[2]
    wi=array(rep(0,M*M*10),c(M,M,10))
    for (i in 1:10){
      GLP = huge(S, lambda = rholist[i], method = "glasso", cov.output = TRUE, verbose = TRUE)
      #print(GLP$icov[[1]])
      wi[,,i] = as.matrix(GLP$icov[[1]])
    }
    #GLP = glassopath((S), rholist ,approx=TRUE, w.init=S , wi.init=solve(S) ,trace=0)
    #print(wi)
    #print(dim(wi))
    nonzeros = apply(wi, 3, function(x) mean(x[upper.tri(x)]!=0))
    #print(wi)
    #print(dim(wi))
    max.rho = max(rholist[nonzeros > 0.05], rholist[1])
    ## Now make the list of rhos
    rholist = seq(low_range, max.rho, length=10)
    #print("glasso")
    #print(rholist)
 #  }
   
  ## to avoid seeing NAN in mu becuase of having all w=0 in one group we need to
  ## get read of all the rows in ts with w = 0
  ## to avoid seeing NAN in mu becuase of having all w=0 in one group we need to
  ## get read of all the rows in ts with w = 0
  ts_new = ts[W!=0,]
  W_new = W[W!=0]
  #print(rholist)
  if(length(W_new)>(2*k)){
    n = nrow(as.matrix(ts_new))
    folds = cvFolds(n, k, type="consecutive")
    loglike = matrix(0, length(rholist) , k)
    for (ki in 1:k) {
      mu = compute_mu_mixture(W_new[folds$which!=ki], ts_new[folds$which!=ki,])
      S_train = compute_sigma_mixture_emp(W_new[folds$which!=ki] ,ts_new[folds$which!=ki,], mu)
      mu = compute_mu_mixture(W_new[folds$which==ki] , ts_new[folds$which==ki,])
      S_test = compute_sigma_mixture_emp(W_new[folds$which==ki] ,ts_new[folds$which==ki,], mu)
        for(i in 1:length(rholist)){
          # 
          GLP = huge((S_train), lambda = rholist[i] , method = "glasso", cov.output = TRUE, verbose = TRUE)
          #GLP = glasso((S_train), rholist[i] ,  start="warm" , w.init=S_train , wi.init=solve(S_train) )
          # 
          loglike[i, ki] = log_likelihood(as.matrix(GLP$icov[[1]]), S_test)
        }
    }
  
  ind = which.max(rowMeans(loglike))
  rhomax = rholist[ind]
  }else{rhomax = low_range}
  #print("rhomax")
  #print(rhomax)
return(rhomax)
}

log_likelihood <- function(precision, emp_cov) {
  p      <- nrow(precision)
  logdet <- determinant(precision, logarithm=T)$modulus
  loglik <- 0.5 * (logdet - sum(emp_cov * precision) - p*log(2*pi))
  return(as.numeric(loglik))
}
compute_mu_mixture<-function(W,lambda){
  N = dim(lambda)[1]
  d = dim(lambda)[2]
  muSum = matrix(0, 1, d)
  WSum = 0
  for (i in 1:N){
    muSum = muSum + W[i]*lambda[i,]
    WSum = WSum + W[i]
  }
return(muSum/WSum) 
}
compute_sigma_mixture_emp<-function(W, lambda, mu){
  N = dim(lambda)[1]
  d = dim(lambda)[2]
  SigmaSum = matrix(0, d, d)
  WSum = 0
  for (i in 1:N){
    SigmaSum = SigmaSum + W[i]*t(lambda[i,] - mu)%*%(lambda[i,] - mu)
    WSum = WSum + W[i]
  }
  Sigma_emp = SigmaSum/WSum 
  Sigma_emp = make.positive.definite(Sigma_emp,tol=1e-5)
 
  return(make.positive.definite(Sigma_emp ,tol=1e-5))
}
 
#################################################################
#' A function to 
#' @param wij
#' @param data  
#' @param K
#' @param invSigma_init
#' @param penalty
#' @return invSigma
#' @export
compute_invSigma_mix_huge<-function(wij,data,K,invSigma_init,penalty){
  N=dim(data)[1]
  M=dim(data)[3]
  invSigma=array(rep(0,K*M*M),c(K,M,M))
  lcluster=apply(wij,1,which.max)
  print("lcluster")
  print(lcluster)
  for(k in 1:K){## each K
    id=which(lcluster==k)
 
    n = as.numeric(table(lcluster)[k])
 
    C=cov.wt(x=data[1:N,k,1:M],wt=wij[1:N,k],method="ML")
    S = make.positive.definite(C$cov,tol=1e-5)
    ##################################### 
    ## STARS
    if(penalty==1){
 
      if (!is.na(n)){
        if(n>M){  #0000000000000000000000000000

          ##############huge
          est = huge(S, nlambda = 10, method = "glasso",lambda.min.ratio=0.005, cov.output = TRUE, verbose = TRUE)
          lam = est$lambda
          print("huge")
          print(lam)
          ##############glasso
low_range = 0.001
rholist = seq(low_range, max(abs(S))+low_range , length=11)[-1]
wi=array(rep(0,M*M*10),c(M,M,10))
for (i in 1:10){
GLP = huge(S, lambda = rholist[i], method = "glasso", cov.output = TRUE, verbose = TRUE)
#print(GLP$icov[[1]])
wi[,,i] = as.matrix(GLP$icov[[1]])
}
#GLP = glassopath((S), rholist ,approx=TRUE, w.init=S , wi.init=solve(S) ,trace=0)
#print(wi)
#print(dim(wi))
nonzeros = apply(wi, 3, function(x) mean(x[upper.tri(x)]!=0))
#print(wi)
#print(dim(wi))
max.rho = max(rholist[nonzeros > 0.05], rholist[1])
## Now make the list of rhos
rholist = seq(low_range, max.rho, length=10)
lam = rholist
          print("glasso")
          print(lam)
          invSigma[k,,] = huge_estimation(lam,  data[,k,1:M]  ) #00000000000000000000000000
          invSigma[k,,] = huge_estimation(lam, t(t(data[,k,1:M])%*%diag(wij[1:N,k]))) #00000000000000000000000000
          #print('new K')
          for (i in 1:M){
            if(is.na(invSigma[k,i,i]) || is.infinite(invSigma[k,i,i]) ){ invSigma[k,i,i]= 0}
          }
         
        }
      }else{ #if(n>M){
        ro = 0.001
        #print(c("ro",ro))
        est = huge(data[id,k,1:M], lambda = ro, method = "glasso", cov.output = TRUE, verbose = TRUE)
        invSigma[k,,] = as.matrix(est$icov[[1]]) 
        #print('new K')
        for (i in 1:M){
          if(is.na(invSigma[k,i,i]) || is.infinite(invSigma[k,i,i])){ invSigma[k,i,i]= 0}
        }
         
      }
    }#penalty 1
    #####################################
    ## find the rho by 2(M)^2/(N||invSigma||)
    if(penalty==2){
      ro = 2*(M^2)/(N*sum(abs(invSigma_init[k,,])))
      #print(c("ro",ro))
      est = huge(data[id,k,1:M], lambda = ro, method = "glasso", cov.output = TRUE, verbose = TRUE)
      invSigma[k,,] = as.matrix(est$icov[[1]])
    }#penalty 2
    #####################################
    ## different rho in each iteration
    if(penalty==3){
      if(n>2){
        ro = 2*(M^2)/(n*sum(abs(solve(make.positive.definite(cov(lambda[id,k,1:M]),tol=1e-5))))) 
        #print(c("ro",ro))
      }else{
        ro = 0.01
      } 
      est = huge(data[id,k,1:M], lambda =ro, method = "glasso", cov.output = TRUE, verbose = TRUE)
      invSigma[k,,] = make.positive.definite(as.matrix(est$icov[[1]]),tol=1e-5)

    }#penalty 3
    #####################################
      for (i in 1:M){
            if(is.na(invSigma[k,i,i]) || is.infinite(invSigma[k,i,i]) ){ invSigma[k,i,i]= 0}
    }
    invSigma[k,,] = make.positive.definite(invSigma[k,,] ,tol=1e-5)
  }#each K
  return(invSigma)
} 
####################################################################
#' A function to    
#'    
#'    
#' @param lam 
#' @param X1
#' @return as.matrix(inv_sigma_huge)
#' @export
huge_estimation <- function(lam, X1){
  #print(X1)
  write.csv(X1, "X.csv")
  est = huge(X1, nlambda = 10, method = "glasso",lambda.min.ratio=0.005, cov.output = TRUE, verbose = TRUE)
  print("huge estimation")
  print(est$lambda)
  #a = which.max(est$loglik)
  final = huge.select(est, criterion = "stars", ebic.gamma = 0.5, stars.thresh = 0.1,
                      stars.subsample.ratio = NULL, rep.num = 20, verbose = TRUE)

  inv_sigma_huge = final$opt.icov
  #inv_sigma_huge = est$icov[[a]]
 
  return(as.matrix(inv_sigma_huge))
}
  
###########################################
#' A function to    
#'   
#' @param CS 
#' @param NS
#' @return val
#' @export
myfr<-function(CS,NS ){
  if(  NROW(CS)!=NCOL(CS) || NROW(NS)!=NCOL(NS) || NROW(CS)!=NROW(NS) ){
    stop ("Error in myfr function")
  }
  Xp=CS
  Yp=NS
  Xp[lower.tri(Xp,diag=FALSE)]=0
  Yp[lower.tri(Yp,diag=FALSE)]=0
  n=(NROW(Xp)*(NROW(Xp)+1))/2
  
  z1=abs(Xp-Yp)/(abs(Xp)+1e-2)
 
  val=sum(z1)/n
   
  return(val)
}

############################################
#' MixMPLN main
#'
#' This function receives the sample-taxa count matrix and cluser the samples based on the MixMPLN method.  
#' @param data
#' @param K number of components
#' @param penalty "no_sparsity", "CV", "StARS" , "fixed" , "iterative"
#' @param init "Kmeans"
#' @param rep number of repeat
#' @param out "precision" , "partial" , "adj"
#' @param threshold a value between 0 and 1 which indicates the threshod to generate the adjacency matrix
#' @return
#' @export


MixMPLN <- function(X,K, penalty , init , rep ){

if( init=="Kmeans"){
  rep = 1
}#if( init=="Kmeans")

out = MixMPLN_function(X,K, penalty , init)
 
ll_print = list()
ll_print[1] = out$ll
if(rep>1){
  for (i in 2:rep){
    out_new = MixMPLN_function(X,K, penalty , init)
    ll_print[i] = out_new$ll
    if(out_new$ll > out$ll){
      out = out_new
    }#if(out_new$ll)
}#for (i in 1:rep)
}#if(rep>1)


return(out)
}# function MixGGm
###############################################
