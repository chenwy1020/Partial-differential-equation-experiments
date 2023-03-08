n=uniform(9,100,1000);
for i = 1:1:10
    m=i*100;
    y(i)=SimpsonCS2(m,0,1);
    h(i)=1/m;
    H(i)=power(h(i),4);
end
for i= 1:1:10
    r(i)=abs(y(i)/(2.0/3)-1);
    R(i)=2.0/3-y(i);
end