function ce= LDA_clime(s0,p)
n1=200;
n2=200;
n=n1+n2;
%p=800;
mu1=zeros(p,1);
mu2=zeros(p,1);
%s0=10;
mu2(1:s0)=1;

%Model 1
% Sigma=ones(p,p)*0.5;
% for i=1:p
%    Sigma(i,i)=1;
% end


%Model 2
 for i=1:p
     for j=1:p
         Sigma(i,j)=0.8^(abs(i-j));
     end
  end


%Model 3
%BB = zeros(p,p);
%for i = 1:p
%  for j = 1:p
%    if i == j
%      BB(i, j) = 1;
%    elseif i <=s0 & i< j
%      BB(i, j) = 0.5*binornd(1,0.2); 
%    elseif (s0+1)<=i & i<j
%      BB(i, j)= 0.5;
%    else
%      BB(i, j) = BB(j,i);
%    end
%  end
% end
%lam_min = min(eig(BB));
%deldel= max(-lam_min,0)+0.05;
%%Rho = (BB+10*deldel*ones(p, p))/(1+deldel);
%Rho = (BB+10*deldel*eye(size(BB)))/(1+deldel);
%Sigma = inv(Rho);


%Model 4
% for i=1:p
%     for j=1:p
%         Rho(i,j)=0.2^(abs(i-j));
%     end
%  end
%Sigma = inv(Rho);

%SIM=10;
%ce=zeros(SIM,1);

%hs=Sigma^0.5;

%for si=1:SIM
%x=mvnrnd(mu1,Sigma,n1);
%x=contaminate_data(x,0.1);
x=mvt(mu1',Sigma,2,n1); %自由度为5的t分布
%x = mixed_sample(n1,0.8,mu1,Sigma); %混合分布
%x=mvt(mu1',Sigma,1,n1); %自由度为1的t分布（柯西分布）

%y=mvnrnd(mu2,Sigma,n2);
%y=contaminate_data(y,0.1);
y=mvt(mu2',Sigma,2,n2); %自由度为5的t分布
%y = mixed_sample(n2,0.8,mu2,Sigma); %混合分布
%y=mvt(mu2',Sigma,1,n2); %自由度为1的t分布（柯西分布）

s1=cov(x);
s2=cov(y);
st=(s1*n1+s2*n2)/n;
smx=mean(x,1);
smy=mean(y,1);
smxy=smx-smy;
tmxy=(smx+smy)/2;




%%cross-validataion
K=10;


x0=x;
y0=y;
LL=10;
cv=zeros(LL,1);

for ei=1:LL
epsilon=0.01*ei;
cv(ei)=0;
for ki=1:K
    x=x0;
    locix=(1+(ki-1)*n1/K):(n1/K*ki);
    %loci=(1+(ki-1)*n1/K):(n1/K*ki);
    %x1=x(loci,:);
    x1=x(locix,:);
    x0=x;
    %x(loci,:)=[];
    x(locix,:)=[];
    x2=x;

    y=y0;
    lociy=(1+(ki-1)*n2/K):(n2/K*ki);
    %y1=y(loci,:);
    y1=y(lociy,:);
    y0=y;
    %y(loci,:)=[];
    y(lociy,:)=[];
    y2=y;
     
    for i=1:(n1/K)
    cv(ei)=cv(ei)+ldac(x0(locix(i),:),x2,y2,epsilon);
    end

    for i=1:(n2/K)
    %cv(ei)=cv(ei)+1-ldac(y0(locix(i),:),x2,y2,epsilon);
    cv(ei)=cv(ei)+1-ldac(y0(lociy(i),:),x2,y2,epsilon);
    end
end

cv(ei)=cv(ei)/(n1+n2);

end
[cvm,locm]=max(cv);

%计算分类误差

% 设置参数：精度和迭代次数
pdtol = 1e-3;
pdmaxiter = 50;
cgtol = 1e-8;
cgmaxiter = 200;
epsilon=0.01*locm;
%x0=inv(st+0.1*eye(p))*smxy;
x0=(st+0.1*eye(p))/smxy;
beta = clime(x0,st, smxy', epsilon, pdtol, pdmaxiter, cgtol, cgmaxiter);


%z=mvnrnd(mu1,Sigma,n1);
%z=contaminate_data(z,0.1);
z=mvt(mu1',Sigma,2,n1); %自由度为3的t分布
%z = mixed_sample(n1,0.8,mu1,Sigma); %混合分布
%z=mvt(mu1',Sigma,1,n1); %自由度为1的t分布（柯西分布）

result1=zeros(n1,1);
for i=1:n1
    result1(i,1)=((z(i,:)-tmxy)*beta>0)*1;
end

%z=mvnrnd(mu2,Sigma,n2);
%z=contaminate_data(z,0.1);
z=mvt(mu2',Sigma,2,n2); %自由度为3的t分布
%z = mixed_sample(n2,0.8,mu2,Sigma); %混合分布
%z=mvt(mu2',Sigma,1,n2); %自由度为1的t分布（柯西分布）


result2=zeros(n2,1);
for i=1:n2
    result2(i,1)=1-((z(i,:)-tmxy)*beta>0)*1;
end

%ce(si)=1-mean(result1+result2)/2;
ce=1-mean(result1+result2)/2;
end

