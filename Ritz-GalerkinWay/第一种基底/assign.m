function [y] = assign(n,x,c,N)

for i = 1:1:(n+1)
    u=0;
    for j=1:1:N
        u = u+c(j)*sin(j*pi*x(i));
    end
    y(i)=power(sin(x(i)/sin(1)-x(i)-u),2);
end
end
