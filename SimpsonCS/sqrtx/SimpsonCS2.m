%% Simpson积分


function csans = SimpsonCS2(n,left,right)
    x = uniform(n,left,right);
    y = assign2(n,x);
    m=n/2;
    temp=0;
    for i = 1:1:m
        temp=temp+(y(2*i-1)+4*y(2*i)+y(2*i+1))*(x(i+1)-x(i))/3;
    end
    csans = temp;
    
end