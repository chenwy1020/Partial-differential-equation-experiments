function [I]=ComGauss(X)

z=zeros(3,1);
y=zeros(1,3);
y=[5,8,5]/9;

for i=1:1:3
    z(i)= X(i, : ) * y';
end

I=y*z;

end