function [I]=ComGauss3(X)

z=zeros(3,1);
y=zeros(3,1);
y(1)=5/9; y(2)=8/9; y(3)=5/9;

for i=1:1:3
    z(i)= X(i, : ) * y;
end

I=y'*z;

end