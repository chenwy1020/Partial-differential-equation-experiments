function I = Gauss(a,b,xi,hi,n)
% 三点Gauss积分公式
% n=0,对f*hi积分
% n=1,对q*(2x-1)(x-1)*hi积分
% n=2,对q*4x(1-x)*hi积分
% n=3,对q*(2x-1)x*hi 积分

w=zeros(3,1);
yi=zeros(3,1);

w=[5 8 5]'/9.0;
yi(1)=(a+b)/2.0-(b-a)/2.0*sqrt(3/5.0);
yi(2)=(a+b)/2.0;
yi(3)=(a+b)/2.0+(b-a)/2.0*sqrt(3/5.0);

q=@(x)(pi^2/4.0);

if n==0
    f=@(x)(sin(pi*(xi+hi*x)/2)*pi^2/2*hi);
end

if n==1
    f=@(x)(q(xi+hi*x)*(2*x-1)*(x-1)*hi);
end

if n==2
    f=@(x)(q(xi+hi*x)*4*x*(1-x)*hi);
end

if n==3
    f=@(x)(q(xi+hi*x)*(2*x-1)*x*hi);
end

I=(b-a)/2.0*[f(yi(1)),f(yi(2)),f(yi(3))]*w;

end

