n=uniform(9,100,1000);
for i = 1:1:10
    m=i*100;
    y(i)=SimpsonCS1(m,0,2*pi);
    h(i)=2*pi/m;
    H(i)=power(h(i),4);
end
for i= 1:1:10
    r(i)=abs(y(i)/(-2*pi)-1);
    R(i)=-2*pi-y(i);
end