CODA<-function(s0,p){
  
  library(pracma)
  library(mvtnorm)
  library(glmnet)
  library(LaplacesDemon)
  
  source('F_tu.R')
  source('f_tranf.R')
  #source('Beta_npn.R')
  source('mvt.R')
  source('MixNormal.R')
  source('cv.dsda.R')
  source('cv.folds.R')
  source('dsda.path.R')
  source('predict.dsda.R')
  
  
  n1=200;
  n2=200;
  n=n1+n2;
  #p=100;
  mu1=zeros(p,1);
  mu2=zeros(p,1);
  #s0=10;
  mu2[1:s0]=1;
  
  
  #Model 1
  #Sigma=matrix(0.5,p,p)
  #for (i in 1:p){
  #  Sigma[i,i]=1
  #}
  
  #Model 2
  Sigma=matrix(ncol=p,nrow=p)
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j]=0.8^(abs(i-j))
    }
  }
  
 

  #x0=rmvnorm(n=n1,mean=mu1,sigma=Sigma)
  #x0 <- mvt(mu1, Sigma, nu=2, n1) #自由度为2的多元t分布
  #x0 <- MixNormal(n1,p,0.8,mu1,Sigma)#混合正态分布
  x0 <- mvt(mu1, Sigma, nu=1, n1) #自由度为1的多元t分布(柯西分布)
  
  #y0=rmvnorm(n=n2,mean=mu2,sigma=Sigma)
  #y0 <- mvt(mu2, Sigma, nu=2, n2) #自由度为5的多元t分布
  #y0 <- MixNormal(n2,p,0.8,mu2,Sigma) #混合正态分布
  y0 <- mvt(mu2, Sigma, nu=1, n2) #自由度为1的多元t分布(柯西分布)
  
 
  # estimators
  smx = apply(x0,2,mean)
  smy = apply(y0,2,mean)
  smxy = smx - smy
  tmxy=(smx+smy)/2
  xy0=rbind(x0,y0)
  smu = apply(xy0,2,mean)
  
  S1_hs = cov(x0,method='kendall')
  S2_hk = cov(y0,method='kendall')
  S_h = (n1*S1_hs+n2*S2_hk)/n
  
  #数据转换
  xy=matrix(nrow = n, ncol = p) 
  
  for (i in 1:dim(xy0)[1]){ 
    xy[i,] = f_tranf(x=xy0[i,],xy=xy0,n1=n1,S=S_h,mx=smx,my=smy)
  }
  
  
  x = xy[1:n1,]
  y = xy[(n1+1):n,]
  Y<- c(rep(0,n1),rep(1,n2))
  
  K=10
  
  la<-cv.dsda(xy,Y,K=K) 
  s<-la$bestlambda
  
  obj.path<-dsda.path(xy,Y,lambda=s)
  beta<-obj.path$beta
  
  
  #测试集
  #z10=rmvnorm(n=n1,mean=mu1,sigma=Sigma)
  #z10 <- mvt(mu1, Sigma, nu=2, n1) #自由度为5的多元t分布
  #z10 <-  MixNormal(n1,p,0.8,mu1,Sigma)
  z10 <- mvt(mu1, Sigma, nu=1, n1) #自由度为1的多元t分布(柯西分布)
  
  #z20=rmvnorm(n=n2,mean=mu2,sigma=Sigma)
  #z20 <- mvt(mu2, Sigma, nu=2, n2) #自由度为5的多元t分布
  #z20 <- MixNormal(n2,p,0.8,mu2,Sigma) #混合正态分布
  z20 <- mvt(mu2, Sigma, nu=1, n2) #自由度为1的多元t分布(柯西分布)
  

  pred1<-as.matrix(cbind(1,z10)%*%beta)
  pred1<-ifelse(pred1<0,1,0)
  

  

  pred2<-as.matrix(cbind(1,z20)%*%beta)
  pred2<-ifelse(pred2>0,1,0)
  
  
  ce=1-mean(pred1+pred2)/2
  
  return(ce)
}


