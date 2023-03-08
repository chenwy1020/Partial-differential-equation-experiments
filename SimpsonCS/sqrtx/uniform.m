function [x] = uniform(n,left,right)
h=(right-left)/n;
for i = 1:1:(n+1)
    x(i)=left+h*(i-1);
end
end