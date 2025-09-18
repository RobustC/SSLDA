LS_LDA<-function(s0,p){

library(pracma)
library(mvtnorm)
library(glmnet)
library(LaplacesDemon)

n1=200;
n2=200;
n=n1+n2;
#p=200;
mu1=zeros(p,1);
mu2=zeros(p,1);
#s0=10;
mu2[1:s0]=1;

source('mvt.R')
source('cv.dsda.R')
source('cv.folds.R')
source('dsda.path.R')
source('predict.dsda.R')
source('MixNormal.R')


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


x=rmvnorm(n=n1,mean=mu1,sigma=Sigma)
#x <- mvt(mu1, Sigma, nu=2, n1) #自由度为2的多元t分布
#x <- MixN(n1,p,mu=mu1,Sigma=Sigma,Kappa=0.8)#混合正态分布
#x <- mvt(mu1, Sigma, nu=1, n1) #多元柯西分布


y=rmvnorm(n=n2,mean=mu2,sigma=Sigma)
#y <- mvt(mu2, Sigma, nu=2, n2) #自由度为5的多元t分布
#y <- MixN(n1,p,mu=mu1,Sigma=Sigma,Kappa=0.8) #混合正态分布
#y <- mvt(mu2, Sigma, nu=1, n2) #自由度为1的多元t分布(柯西分布)

# estimators
xy = rbind(x,y)
Y<- c(rep(0,n1),rep(1,n2))

K=10

la<-cv.dsda(xy,Y,K=K) 
s<-la$bestlambda

obj.path<-dsda.path(xy,Y,lambda=s)
beta<-obj.path$beta


z1=rmvnorm(n=n1,mean=mu1,sigma=Sigma)
#z1 <- mvt(mu1, Sigma, nu=2, n1) #自由度为5的多元t分布
#z1 <-  MixN(n1,p,mu=mu1,Sigma=Sigma,Kappa=0.8)
#z1 <- mvt(mu1, Sigma, nu=1, n1) #自由度为1的多元t分布(柯西分布)

pred1<-as.matrix(cbind(1,z1)%*%beta)
pred1<-ifelse(pred1<0,1,0)

#pred1=t(c(1,z1[1,]))%*%beta

z2=rmvnorm(n=n2,mean=mu2,sigma=Sigma)
#z2 <- mvt(mu2, Sigma, nu=2, n2) #自由度为2的多元t分布
#z2 <- MixN(n1,p,mu=mu1,Sigma=Sigma,Kappa=0.8) #混合正态分布
#z2 <- mvt(mu2, Sigma, nu=1, n2) #自由度为1的多元t分布(柯西分布)

pred2<-as.matrix(cbind(1,z2)%*%beta)
pred2<-ifelse(pred2>0,1,0)


ce=1-mean(pred1+pred2)/2
#ce = 1-(n1*mean(pred1)+mean(pred2)*n2)/n

return(ce)
}



