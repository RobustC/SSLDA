function Mixed_sample = mixed_sample(n,Kappa,mu,S) 
  
  %混合分布的参数
  weights = [Kappa, 1-Kappa];  % 每个分布的权重
  %means = [0, 5];  % 正态分布的均值

  % 生成样本
  n0 = round(weights(1) * n);
  sample_1 = mvnrnd(mu,S,n0)/(sqrt(Kappa+9*(1-Kappa)));
  sample_2 = mvnrnd(mu,9*S,n-n0)/(sqrt(Kappa+9*(1-Kappa)));
  
  % 合并两个样本
  Mixed_sample = [sample_1; sample_2];
  

  % 打乱样本顺序
  %Mixed_sample = Mixed_sample(randperm(length(Mixed_sample)),:);
  Mixed_sample = Mixed_sample(randperm(n),:);
  
  
end

