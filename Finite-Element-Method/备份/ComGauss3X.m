function [I]=ComGauss3X(X)

y=zeros(3,1);
y(1)=5/9; y(2)=8/9; y(3)=5/9;
I= y'*X;

end