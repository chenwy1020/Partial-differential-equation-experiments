n=uniform(99,100,10000);
for i = 1:1:100
    m=i*100;
    y(i)=SimpsonCS(m,0,2*pi);
    h(i)=1/m;
    H(i)=power(h(i),4);
end
for i= 1:1:100
    r(i)=abs(y(i)/(-2*pi)-1);
    R(i)=-2*pi-y(i);
end
xlswrite('D:\Documents\MATLAB\SimpsonCS\date_xsinx_H_R.xlsx',H,'A1:A100')
xlswrite('D:\Documents\MATLAB\SimpsonCS\date_xsinx_H_R.xlsx',R,'B1:B100')