%切片积分
function [I]=ComSimpson(X)

[N,~]=size(X);
N=N-1;
h=1/N;
y=2*ones(1,N+1);

y(1, 1)=1; y(1, 2:2:N)=4; y(1, N+1)=1;

I=y*X*y'*h*h/3/3;

end

