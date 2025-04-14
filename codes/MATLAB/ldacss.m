function czss=ldacss(z,x,y,epsilon)
n1=size(x,1);
n2=size(y,1);
n=n1+n2;
p=size(x,2);
smx=spatialMedian(x);
smy=spatialMedian(y);
s1=spatialSignCovariance(x,smx);
s2=spatialSignCovariance(y,smy);
st=(s1*n1+s2*n2)/n;
smxy=smx-smy;
tmxy=(smx+smy)/2;

normDiff=0;
    for i = 1:n1
        diff = x(i, :) - smx; % 数据点相对于中位数的差值
        normDiff =normDiff+norm(diff);   % 差值的范数
    end

    for i = 1:n2
        diff = y(i, :) - smy; % 数据点相对于中位数的差值
        normDiff =normDiff+norm(diff);   % 差值的范数
    end
normDiff=normDiff/n;
st=st*normDiff^2;
% 设置参数：精度和迭代次数
pdtol = 1e-3;
pdmaxiter = 50;
cgtol = 1e-8;
cgmaxiter = 200;
%epsilon=0.01;
%x0=inv(st+0.1*eye(p))*smxy;
x0=(st+0.1*eye(p))/smxy;
beta = clime(x0,st, smxy', epsilon, pdtol, pdmaxiter, cgtol, cgmaxiter);
czss=1*((z-tmxy)*beta>0);
end