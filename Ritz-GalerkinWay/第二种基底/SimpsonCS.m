function csans = SimpsonCS(n,x,y)
    temp=0;
    for i = 1:1:(n-3)
        temp=temp+(y(i+2)+4*y(i+1)+y(i))*(x(i+1)-x(i))/6;
    end
    csans = temp;
    
end


