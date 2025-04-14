function cz=ldac(z,x,y,epsilon)
n1=size(x,1);
n2=size(y,1);
n=n1+n2;
p=size(x,2);
s1=cov(x);
s2=cov(y);
st=(s1*n1+s2*n2)/n;
smx=mean(x,1);
smy=mean(y,1);
smxy=smx-smy;
tmxy=(smx+smy)/2;
% 设置参数：精度和迭代次数
pdtol = 1e-3;
pdmaxiter = 50;
cgtol = 1e-8;
cgmaxiter = 200;
%epsilon=0.01;
%x0=inv(st+0.1*eye(p))*smxy;
x0=(st+0.1*eye(p))/smxy;
beta = clime(x0,st, smxy', epsilon, pdtol, pdmaxiter, cgtol, cgmaxiter);
cz=1*((z-tmxy)*beta>0);
end