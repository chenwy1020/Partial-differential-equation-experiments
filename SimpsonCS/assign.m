function [y] = assign(n,x)

for i = 1:1:(n+1)
    y(i) = x(i)*sin(x(i));
end
end
