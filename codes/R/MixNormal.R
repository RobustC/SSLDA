MixNormal <- function(n,p,Kappa,mu,S) {
  
  #混合分布的参数
  weights = c(Kappa, 1-Kappa)  # 每个分布的权重
  #means = rep(0,p)  # 正态分布的均值
  
  # 生成混合分布数据
  sample_size_1 <- rmvnorm(n=round(weights[1] * n),mean=mu,sigma=S)/(sqrt(Kappa+9*(1-Kappa)))
  sample_size_2 <- rmvnorm(n=n-round(weights[1] * n),mean=mu,sigma=9*S)/(sqrt(Kappa+9*(1-Kappa)))
  
  # 合并两个样本
  mixed_sample <- rbind(sample_size_1, sample_size_2)
  
  # 打乱顺序
  mixed_sample <- mixed_sample[sample(nrow(mixed_sample)), ]
  
  return(mixed_sample)
  
}