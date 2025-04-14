mvt=function(mu, Sigma, nu, n){
  # 生成多元 t 分布数据
  # mu: 均值向量 (1 x d)
  # Sigma: 协方差矩阵 (d x d)
  # nu: 自由度
  # n: 样本数量
  # 返回: n x d 的多元 t 分布样本
  
  d = length(mu) # 数据维度
  
  # 生成标准多元正态分布样本
  Z = rmvnorm(n=n,mean=zeros(1, d),sigma=Sigma)
  
  # 生成卡方分布样本，并归一化
  chi2Samples = rchisq(n=n, df=nu)
  scalingFactors = sqrt(nu / chi2Samples)
  
  # 生成多元 t 分布样本
  T <- t(apply(Z * scalingFactors, 1, function(x) x + mu))
  
}

