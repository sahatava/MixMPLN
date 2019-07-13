
library(ROCR)


##########################################################
#' A function to    
#'    
#' @param real
#' @param res
#' @return val
#' @export
frob <- function(CS,NS){
  if( NROW(CS)!=NCOL(CS) || NROW(NS)!=NCOL(NS) || NROW(CS)!=NROW(NS) ){
    stop ("Error in myFrobenius function")
  }
  Xp = -cov2cor(CS)
  Yp = -cov2cor(NS)
  z = Xp-Yp
  val = sqrt(sum(z*z)) 
 
  return(val)
}
###########################################
#' A function to    
#'    
#' @param real
#' @param res
#' @return val
#' @export
myfr2<-function(CS,NS ){
  if(  NROW(CS)!=NCOL(CS) || NROW(NS)!=NCOL(NS) || NROW(CS)!=NROW(NS) ){
    stop ("Error in myfr function")
  }
  Xp=CS
  Yp=NS
  Xp[lower.tri(Xp,diag=FALSE)]=0
  Yp[lower.tri(Yp,diag=FALSE)]=0
  n=(NROW(Xp)*(NROW(Xp)+1))/2
  
  z1= (abs(Xp-Yp)/(abs(Xp)+1e-2) + abs(Xp-Yp)/(abs(Yp)+1e-2) )/2
 
  val=sum(z1)/n
   
  return(val)
}
############################################
#' A function to    
#'    
#' @param real
#' @param res
#' @return val
#' @export
myrms<-function(CS,NS){
  if( NROW(CS)!=NCOL(CS) || NROW(NS)!=NCOL(NS) || NROW(CS)!=NROW(NS) ){
    stop ("Error in myrms function")
  }
  Xp=CS
  Yp=NS
  Xp[lower.tri(Xp,diag=FALSE)]=0
  Yp[lower.tri(Yp,diag=FALSE)]=0
  n=(NROW(Xp)*(NROW(Xp)+1))/2
   
  z1=Xp-Yp
   
  val=sqrt(((sum(z1*z1)) )/n)
#  print(c("val",val))
  return(val)
}
#########################################
#' A function to    
#'    
#' @param real
#' @param res
#' @return val
#' @export
mysnsp<-function(inv,res){
  

  thresh = 0
  if(NROW(inv)!=NCOL(inv) || NROW(res)!=NCOL(res) || NROW(inv)!=NROW(res)){
    stop ("Error in mysnsp function")
  }
  Tc=(abs(cov2cor(as.matrix(inv))) >= thresh)
  Pc=(abs(cov2cor(as.matrix(res))) >= thresh)
  n=NROW(inv)
  TP=0
  FP=0
  FN=0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(Tc[i,j]==TRUE && Pc[i,j]==TRUE){
        TP=TP+1
      }else{
        if(Tc[i,j]==TRUE && Pc[i,j]==FALSE){
          FN=FN+1
        }else{
          if(Tc[i,j]==FALSE && Pc[i,j]==TRUE){
            FP=FP+1
          }
        }
      }
    }
  }
  if(TP!=0){
    sn=TP/(TP+FN)
    sp=TP/(TP+FP)
    fv=(2*sn*sp)/(sn+sp)
  }else{
    sn = sp = fv = 0
  }
  
  return(list("sn"=sn,"sp"=sp,"fv"=fv))
}
#########################################################
#' A function to    
#'    
#' @param real
#' @param res
#' @return res
#' @export
cor2directed <- function(input, ts){
input = cov2cor(input)
res = input
for(i in 1:nrow(input)){
  for(j in 1:ncol(input)){
    if(input[i,j]> ts){res[i,j]=1}
    if(input[i,j]< -ts){res[i,j]=1}
    if(input[i,j]>= -ts && input[i,j]<= ts ){res[i,j]=0}
  }
}
return(res)
}
########################################
#' A function to    
#'    
#' @param real
#' @param res
#' @return val
#' @export
myROC <- function(real, res){

inv1 = cor2directed(real, 0)
labelsA = inv1[upper.tri(inv1, diag=FALSE)]

predictions = abs(res[upper.tri(res, diag=FALSE)])
pred <- prediction(predictions , labelsA  )
area = performance(pred, measure = "auc")
val = slot(area,"y.values")[[1]]

return(val)
}
#######################################

######################################
#' A function to    
#'    
#' @param inv1
#' @param inv2
#' @param res1
#' @param res2
#' @return val
#' @export
data2_split <- function(inv1,inv2,res1,res2){

 
res1 = res[1:d,]
res2 = res[(d+1):(2*d),]
 
index = which.min(c( frob(inv1,res1) , frob(inv1,res2)))
if( index==1 ){
  A = res1
  B = res2
}
if( index==2 ){
  A = res2
  B = res1
}
 
 
 
return(list("A"=A  , "B"=B))
}
################################################################
#' A function to    
#'    
#' @param inv1
#' @param inv2
#' @param inv3
#' @param res1
#' @param res2
#' @param res3
#' @return val
#' @export
data3_split <- function(inv1,inv2,inv3,res1,res2,res3){


res1 = res[1:d,]
res2 = res[(d+1):(2*d),]
res3 = res[((2*d)+1):(3*d),]

index = which.min(c(frob(inv1,res1),frob(inv1,res2),frob(inv1,res3)))
if( index == 1 ){
  A = res1
  if(frob(inv2,res2) < frob(inv2,res3)){
    B = res2
    C = res3
  }else{
    B = res3
    C = res2
  } 
}
if( index == 2 ){
  A = res2
  if(frob(inv2,res1) < frob(inv2,res3)){
    B = res1
    C = res3
  }else{
    B = res3
    C = res1
  }
  
}
if( index == 3 ){
  A = res3
  if(frob(inv2,res1) < frob(inv2,res2)){
    B = res1
    C = res2
  }else{
    B = res2
    C = res1
  }  
}

  
 
return(list("A"=A , "B"=B , "C"=C))
}
##############################
#' A function to    
#'    
#' @param real
#' @param res
#' @return val
#' @export
compare <- function(real , res){

K = dim(as.matrix(res))[1]

if(K == 1){
out1 = frob(real,res)
out2 = myfr2(real,res)
out3 = myrms(real,res)
out4 = mysnsp(real,res)$sn
out5 = mysnsp(real,res)$sp
out6 = myROC(real,res)
}

if(K == 2){
res1 = data2_split(real[[1]],real[[2]],res[[1]],res[[2]])$A
res2 = data2_split(real[[1]],real[[2]],res[[1]],res[[2]])$B

out1 = mean(frob(real[[1]],res1) , frob(real[[2]],res2)  )
out2 = mean(myfr2(real[[1]],res1) , myfr2(real[[2]],res2) )
out3 = mean(myrms(real[[1]],res1) , myrms(real[[2]],res2) )
out4 = mean(mysnsp(real[[1]],res1)$sn , mysnsp(real[[2]],res2)$sn )
out5 = mean(mysnsp(real[[1]],res1)$sp , mysnsp(real[[2]],res2)$sp )
out6 = mean(myROC(real[[1]],res1) , myROC(real[[2]],res2) )
}

if(K == 3){
res1 = data3_split(real[[1]],real[[2]],real[[3]],res[[1]],res[[2]],res[[3]])$A
res2 = data3_split(real[[1]],real[[2]],real[[3]],res[[1]],res[[2]],res[[3]])$B
res3 = data3_split(real[[1]],real[[2]],real[[3]],res[[1]],res[[2]],res[[3]])$C

out1 = mean(frob(real[[1]],res1) , frob(real[[2]],res2) , frob(real[[3]],res3) )
out2 = mean(myfr2(real[[1]],res1) , myfr2(real[[2]],res2) , myfr2(real[[3]],res3) )
out3 = mean(myrms(real[[1]],res1) , myrms(real[[2]],res2) , myrms(real[[3]],res3))
out4 = mean(mysnsp(real[[1]],res1)$sn , mysnsp(real[[2]],res2)$sn , mysnsp(real[[3]],res3)$sn)
out5 = mean(mysnsp(real[[1]],res1)$sp , mysnsp(real[[2]],res2)$sp , mysnsp(real[[3]],res3)$sp)
out6 = mean(myROC(real[[1]],res1) , myROC(real[[2]],res2) , myROC(real[[3]],res3))
}

return( list("frob"=out1 , "fr"=out2 , "rms"=out3 , "sn"=out4 , "sp"=out5 , "ROC"=out6) )
}

 
