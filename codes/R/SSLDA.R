SSLDA_NEW<-function(s0,p,threshold=0.1){
  library(CVXR)
  library(Rmosek)
  library(pracma)
  library(mvtnorm)
  
  library(ICSNP)
  library(SpatialNP)
  
  
  
  source('mvt.R')
  source('MixN.R')

  
  ## 1. 生成数据
  n1=200;
  n2=200;
  n=n1+n2;
  mu1=pracma::zeros(p,1);
  mu2=pracma::zeros(p,1);
  mu2[1:s0]=1;
  
  
  #Model 1
  #Sigma=matrix(0.5,p,p)
  #diag(Sigma)<-1
  
  
  #Model 2
  Sigma=matrix(ncol=p,nrow=p)
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j]=0.8^(abs(i-j))
    }
  }
  
  
  
  #x=rmvnorm(n=n1,mean=mu1,sigma=Sigma)
  #x <- mvt(mu1, Sigma, nu=2, n1) #自由度为2的多元t分布
  #x <- MixN(n1,p,mu=mu1,Sigma=Sigma,Kappa=0.8)#混合正态分布
  x <- mvt(mu1, Sigma, nu=1, n1) #多元柯西分布
  
  
  #y=rmvnorm(n=n2,mean=mu2,sigma=Sigma)
  #y <- mvt(mu2, Sigma, nu=2, n2) #自由度为2的多元t分布
  #y <- MixN(n2,p,mu=mu2,Sigma=Sigma,Kappa=0.8) #混合正态分布
  y <- mvt(mu2, Sigma, nu=1, n2) #自由度为1的多元t分布(柯西分布)
  
  
  
  
  # 空间中位数
  smx = spatial.median(x)
  smy = spatial.median(y)
  tmxy=(smx+smy)/2
  
  

  s1 = SCov(x)
  s2 = SCov(y)
  
  S_h=(s1*n1+s2*n2)/n
  delta_h = smx-smy
  
  
  #true parameters
  #truedalta = mu1-mu2
  #omega = psych::tr(Sigma)/p
  #Gamma = Sigma/omega
  
  #truegammA <- solve(Gamma) %*% truedalta
  #truea<-c(1,2,rep(0,p-2))
  
  
  lambda_val <-threshold 
  
  # 变量
  gammA <- Variable(ncol(Sigma))
  
  # 目标函数：L1范数最小化
  objective <- Minimize(norm1(gammA))
  
  # 约束：最大范数小于 lambda_val
  constraints <- list(norm_inf(p*S_h %*% gammA - delta_h) <= lambda_val)
  
  # 建立问题
  problem <- Problem(objective, constraints)
  
  # 求解
  result <- solve(problem, solver = "ECOS")
  
  # 输出结果
  result$status  # 应该是 "optimal"
  result$value   # 最小值
  gammA_h = result$getValue(gammA)  # 得到 gammA 的估计值
  
  #z1=rmvnorm(n=n1,mean=mu1,sigma=Sigma)
  #z1 = mvt(mu1, Sigma, nu=2, n1) #自由度为2的多元t分布
  #z1 <-  MixN(n1,p,mu=mu1,Sigma=Sigma,Kappa=0.8)#混合正态分布
  z1 <- mvt(mu1, Sigma, nu=1, n1) #自由度为1的多元t分布(柯西分布)
  
  result1=rep(0,n1)
  
  for (i in 1:n1){
    result1[i]=(t(as.matrix((z1[i,]-tmxy)))%*%gammA_h>0)*1;
  }
  
  
  #z2=rmvnorm(n=n2,mean=mu2,sigma=Sigma)
  #z2=mvt(mu2, Sigma, nu=2, n2) #自由度为2的多元t分布
  #z2 <- MixN(n2,p,mu=mu2,Sigma=Sigma,Kappa=0.8) #混合正态分布
  z2 <- mvt(mu2, Sigma, nu=1, n2) #自由度为1的多元t分布(柯西分布)
  
  result2=rep(0,n2)
  
  for (i in 1:n2){
    result2[i]=1-(t(as.matrix((z2[i,]-tmxy)))%*%gammA_h>0)*1;
  }
  
  
  
  ce=1-mean(result1+result2)/2
  
  return(ce)
}