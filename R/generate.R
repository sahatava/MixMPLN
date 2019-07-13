library(mvtnorm)
library(corpcor)
 
 
 
####################################################################
#' A function to  
#'
#' @param nrow 
#' @param sparse
#' @param min 
#' @param max
#' @return M
#' @export
generate_random_pd_sym_matrix<-function(nrow,sparse,min,max){
  if(sparse<0 || sparse>1 || nrow<=0){
    stop ("Error: Sparse needs to be between 0 and 1, and nrow needs to be positive integer")
  }
  M=matrix(0,nrow=nrow,ncol=nrow)
  x=runif(nrow*(nrow-1)/2, min=min,max=max)
  y=sample.int(nrow*(nrow-1)/2, size=sparse*nrow*(nrow-1)/2)
  x[y]=0
  k=1
  for(i in 1:(nrow-1)){
    for(j in (i+1):nrow){
      if(abs(x[k])>1e-2){
        M[i,j]=x[k]
      }else{
        M[i,j]=0
      }
      k=k+1
    }
  }
  M=M+t(M)
  M=M+diag(runif(nrow,min=0.5,max=1.5),nrow=nrow)
  cond=rank.condition(M,tol=1e-5)$condition
  while(cond>30 || !is.positive.definite(M,tol=1e-5)){
    M=M+diag(runif(nrow,min=0.05,max=0.5),nrow=nrow)
    cond=rank.condition(M,tol=1e-5)$condition
  }
  return(M)
}
##########################################################
#' A function to  
#' 
#' @param  X_samp
#' @return list(TMM=X_TMM)
#' @export
normalized_data <- function(X_samp){
  N = nrow(X_samp)
  t <-as.vector(calcNormFactors(as.matrix(t(X_samp)), method = "TMM"))
  X_TMM = X_samp
  for(i in 1:N){
    X_TMM[i,] = X_samp[i,]/t[i]
  }
 
  return(list(TMM=X_TMM))
}
##########################################################
#' A function to    
#' @param y 
#' @param MinValue
#' @param MaxValue 
#' @return x
#' @export
sampling <-  function(y, MinValue , MaxValue){
  N = nrow(y)
  d = ncol(y)
  x=matrix(NA,ncol=d,nrow=N) 
  for (i in 1:N){
    x[i,] = rmultinom(1, size = runif(1,MinValue , MaxValue), prob = y[i,])
  }
	 
  return(x)
}

############ generate data####################

#' A function to generate the synthetic data
#'
#' @param K number of componenet in the synthetic data
#' @param N number of samples in the synthetic sample-taxa count matrix
#' @param d number of taxa in the synthetic sample-taxa count matrix
#' @param sp sparsity level. value between 0 and 1. 
#' @param type "orig" original count data, "samp" sampled data , "TMM" after TMM normalization
#' @return
#' @export
generate <- function(K , N , d , sp, type){
     
     real_precision = list()
   
     invcovmat=generate_random_pd_sym_matrix(nrow=d,sparse=sp,min=-1,max=1)
     #print(invcovmat)
     real_precision[[1]] = invcovmat
     covmat=solve(invcovmat)
     cmean = rep(4,d)
     lamb=rmvnorm(N,mean=cmean,sigma=covmat)
     inv = invcovmat
     
     X = t(apply(lamb,1,function(x) rpois(length(x),exp(x))))
     if(type == "samp") {
       X = sampling(X, 5000 , 10000)
     } 
     if(type == "TMM") {
     samp = sampling(X, 5000 , 10000)
     X = normalized_data(samp)$TMM
     }
   if(K == 1) {
     X = X
   }else {
     for(l in 2:K){
       invcovmat=generate_random_pd_sym_matrix(nrow=d,sparse=sp,min=-1,max=1)
       #print(invcovmat)
       real_precision[[l]] = invcovmat
       covmat=solve(invcovmat)
       cmean = rep(l+2,d)
       lamb=rmvnorm(N,mean=cmean,sigma=covmat)
       inv = invcovmat
       X_new = t(apply(lamb,1,function(x) rpois(length(x),exp(x))))
       if(type == "samp"){
         X_new = sampling(X_new, 5000 , 10000)
       }# end if(type == "samp") 
       if(type == "TMM"){
         samp = sampling(X_new, 5000 , 10000)
         X_new = normalized_data(samp)$TMM
       }# end if(type == "TMM")
     X =rbind(X , X_new)
     }# end for
   }# end else
 
 
  return( list("M"=X, "real_precision"=real_precision))
}

