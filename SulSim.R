rm(list=ls())
cat("\f")

#==============================
# Global settings
#==============================
n=50#10,15,20,30,50,100
p=4
q=4
det="trend"#trend, const
Mu_X=25.027
Sig_X=9.167
max_pq=max(p,q)
nSim=1000
spreadScale=0.1

##############
X=c(28.716265, 26.375666, 30.674826, 28.000638,
    26.290853, 28.346808, 25.444804, 25.304158, 14.697552,
    16.306727, 16.511773, 18.954200, 6.980143, 15.435991,
    13.350322, 18.642814, 15.324487, 16.398606, 24.549812,
    21.063050, 21.431309, 23.936898, 31.115046, 30.086581,
    20.458173, 25.611431, 27.166525, 34.330866, 18.168104,
    9.914513, 10.358615, 24.758457, 17.311942, 22.605556,
    31.593736, 37.570230, 34.901964, 40.737365, 47.157174,
    38.900934, 27.236213, 27.072127, 41.539847, 39.833079)
Y=c(11.928087, 1.476802, 2.136739, -8.647260,
    13.254903, 13.780045, 11.573110, -16.777661, 3.460424,
    -0.642614, -0.737602, -6.007471, -4.707486, -3.204968,
    -1.519991, -3.225722, 1.828321, 1.969864, 0.248302,
    -2.200647, 8.388422, 6.491370, 4.384894, 6.431959,
    2.028613, 5.216008, -0.457179, -10.943970, 0.790508,
    -4.062199, -10.021255, 0.410192, -0.131747, 1.928128,
    8.142723, -3.707299, 11.408643, 15.036579, 17.598418,
    0.694215, 4.568021, 1.939697, 11.083735, -24.345823)

#########
addDeterministic=function(Z,det){
  n=nrow(Z)
  if(det=="none")return(Z)
  if(det=="const")return(cbind(1,Z))
  if(det=="trend")return(cbind(1,1:n,Z))
  stop("det must be none,const,or trend")
}
buildZ=function(Y,X){
  Z=NULL
  for(t in (max_pq+1):length(Y)){
    Zt=NULL
    for(i in 1:p){
      Zt=c(Zt,Y[t-i])
    }
    for(j in 0:q){
      Zt=c(Zt,X[t-j])
    }
    Z=rbind(Z,Zt)
  }
  return(Z)
}
Z=buildZ(Y,X)
Z=addDeterministic(Z,det)
Yc=Y[(max(p,q)+1):length(Y)]
Par=round(coef(lm(Yc~Z-1)),3)
if(det=="const"){
  B0=Par[1]
  Beta_Y=Par[2:(p+1)]
  Beta_X=Par[(p+2):(p+q+2)]
  beta_true=c(B0,Beta_Y,Beta_X)
  k=1+p+(q+1)
  ParNames=c("Intercept",paste0("Y_L",1:p),paste0("X_L",0:q))
}else if(det=="trend"){
  B0=Par[1]
  B_Trend=Par[2]
  Beta_Y=Par[3:(p+2)]
  Beta_X=Par[(p+3):(p+q+3)]
  beta_true=c(B0,B_Trend,Beta_Y,Beta_X)
  k=2+p+(q+1)
  ParNames=c("Intercept","Trend",paste0("Y_L",1:p),paste0("X_L",0:q))
}




Alpha_T=matrix(NA,nSim,k)
Alpha_Q=matrix(NA,nSim,k)

Spred_T=matrix(NA,nSim,k)
Spred_Q=matrix(NA,nSim,k)

library(lpSolve)
library(quadprog)

#==============================
# Build Z matrix
#==============================
buildZ=function(Y,X){
  Z=NULL
  for(t in (max_pq+1):length(Y)){
    Zt=NULL
    for(i in 1:p){
      Zt=c(Zt,Y[t-i])
    }
    for(j in 0:q){
      Zt=c(Zt,X[t-j])
    }
    Z=rbind(Z,Zt)
  }
  return(Z)
}

#==============================
# Data generation (no bounds)
#==============================
simulateFuzzyARDL=function(){
  X=rnorm(n+max_pq,Mu_X,Sig_X)
  Y=numeric(n+max_pq)
  
  for(i in 1:max_pq){
    Y[i]=rnorm(1,Mu_X,Sig_X)
  }
  
  for(t in (max_pq+1):(n+max_pq)){
    Y[t]=B0+
      sum(Beta_Y*Y[(t-1):(t-p)])+
      sum(Beta_X*X[t:(t-q)])+
      rnorm(1,0,1)
  }
  
  Z=buildZ(Y,X)
  Z=addDeterministic(Z,det)
  Yc=Y[(max_pq+1):(n+max_pq)]
  
  Sc=spreadScale*abs(Yc)
  
  return(list(
    Z=Z,
    absZ=abs(Z),
    Yc=Yc,
    Sc=Sc
  ))
}

#==============================
# Tanaka LP estimation
#==============================
fitTanakaLP=function(Z,absZ,Yc,Sc){
  T=nrow(Z)
  k=ncol(Z)
  
  obj=c(rep(0,k),rep(1,k),colSums(absZ))
  
  A1=cbind(Z,-absZ)
  A2=cbind(-Z,-absZ)
  A3=cbind(matrix(0,k,k),-diag(k))
  
  A=rbind(A1,A2,A3)
  rhs=c(Yc+Sc,-(Yc-Sc),rep(0,k))
  dir=rep("<=",nrow(A))
  
  sol=lp("min",obj,A,dir,rhs)
  v=sol$solution
  
  a=v[1:k]
  cpar=pmax(v[(k+1):(2*k)],0)
  
  return(list(a=a,c=cpar))
}

#==============================
# Quadratic Programming estimation
#==============================
makeFullRank=function(Z,tol=1e-10){
  qrZ=qr(Z,tol=tol)
  keep=qrZ$pivot[1:qrZ$rank]
  Z2=Z[,keep,drop=FALSE]
  return(list(Z=Z2,keep=keep,rank=qrZ$rank))
}
fitQP=function(Z,absZ,Yc,Sc,eps=1e-6,tol=1e-10){
  
  fr=makeFullRank(Z,tol=tol)
  
  Zr=fr$Z
  absZr=absZ[,fr$keep,drop=FALSE]
  
  k=ncol(Zr)
  
  Daa=2*crossprod(Zr)
  Dcc=2*diag(k)
  
  Dmat=rbind(
    cbind(Daa,matrix(0,k,k)),
    cbind(matrix(0,k,k),Dcc)
  ) + eps*diag(2*k)
  
  dvec=c(2*crossprod(Zr,Yc),rep(0,k))
  
  G1=t(cbind( Zr, absZr))
  G2=t(cbind(-Zr, absZr))
  G3=t(cbind(matrix(0,k,k),diag(k)))
  
  Amat=cbind(G1,G2,G3)
  bvec=c(Yc+Sc, -(Yc-Sc), rep(0,k))
  
  qp=quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=0)
  
  v=qp$solution
  
  a_r=v[1:k]
  c_r=pmax(v[(k+1):(2*k)],0)
  
  a_full=rep(0,ncol(Z)); a_full[fr$keep]=a_r
  c_full=rep(0,ncol(Z)); c_full[fr$keep]=c_r
  
  return(list(a=a_full,c=c_full,keep=fr$keep,rank=fr$rank,value=qp$value))
}
#==============================
# Performance metrics
#==============================
computeMetrics=function(Yc,YC,SC){
  RMSE=sqrt(mean((Yc-YC)^2))
  MAPE=mean(abs((Yc-YC)/Yc))
  FD=mean(2*SC)
  return(c(RMSE=RMSE,MAPE=MAPE,FD=FD))
}

#==============================
# Monte Carlo simulation
#==============================
Res_T=matrix(NA,nSim,3)
Res_Q=matrix(NA,nSim,3)

for(s in 1:nSim){
  dat=simulateFuzzyARDL()
  
  Z=dat$Z
  absZ=dat$absZ
  Yc=dat$Yc
  Sc=dat$Sc
  
  fitT=fitTanakaLP(Z,absZ,Yc,Sc)
  YC_T=Z%*%fitT$a
  SC_T=absZ%*%fitT$c
  
  fitQ=fitQP(Z,absZ,Yc,Sc)
  YC_Q=Z%*%fitQ$a
  SC_Q=absZ%*%fitQ$c
  
  Alpha_T[s,]=fitT$a
  Alpha_Q[s,]=fitQ$a
  
  Spred_T[s,]=fitT$c
  Spred_Q[s,]=fitQ$c
  
  Res_T[s,]=computeMetrics(Yc,YC_T,SC_T)
  Res_Q[s,]=computeMetrics(Yc,YC_Q,SC_Q)
}

#==============================
# Final results table
#==============================

Res=matrix(NA,k+3,5)
Res[1:k,1]=beta_true
Res[1:k,2]=colMeans(Alpha_T)
Res[1:k,3]=colMeans(Alpha_Q)
Res[1:k,4]=colMeans(Spred_T)
Res[1:k,5]=colMeans(Spred_Q)
Res[(k+1):(k+3),2]=colMeans(Res_T)/sqrt(n)
Res[(k+1):(k+3),3]=colMeans(Res_Q)#/sqrt(n)

rownames(Res)=c(ParNames,"RMSE","MAPE","FD")
colnames(Res)=c("True","Tanaka-a","QP-a","Tanaka-c","QP-c")

Res
write.excel <- function(x,row.names=FALSE,col.names=FALSE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

write.excel(Res[,2:5])
