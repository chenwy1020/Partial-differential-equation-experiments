%仿照老师所给的程序代码，建立Lagrange型二次有限元方法
%此程序为Lagrange型一维二次有限元有限元方法的通用程序
clear

% 进行循环的次数
MK=20;

err2=zeros(MK, 1);      errf=zeros(MK, 1);      
errH1=zeros(MK, 1);   errL2=zeros(MK,1);
CondN=zeros(MK, 1);


for ni=1:1:MK
    n=10*ni;
    a=0; b=1;
    A=zeros(2*n+2, 2*n+2); y=zeros(2*n+2, 1);
    p=zeros(1,2*n+2); I=zeros(4, n);
    
    %% 网格信息
    % 这里考虑均匀剖分 
    h=(b-a)/n;
    p(1,[1:2:2*n+1])=a+[0:1:n]*h;   p(1,[2:2:2*n+2])=a+[0:1:n]*h;
    I(1, : )=[1:2:2*n-1]; I(2, : )=[2:2:2*n];  I(3, : )=[3:2:2*n+1]; I(4, : )=[4:2:2*n+2];

    %% 给出变分问题 a(u, v)
    % 对于此问题 p(x)=1, q(x)=pi^2/4
    p1=1; p2=pi^2/4;

    %% 定义标准区域上的数值求积公式，定义求积节点
    % 这里利用复化Guass公式
    xi(1)=(a+b)/2-(b-a)/2*sqrt(3/5);
    xi(2)=(a+b)/2;
    xi(3)=(a+b)/2+(b-a)/2*sqrt(3/5);
    
    %% 计算I_{i)单元刚度矩阵
    for i=1:1:n
        i1=I(1, i); i2=I(2, i); i3=I(3,i); i4=I(4,i);
        x1=p(i1);  x2=p(i3);  
        hi=x2-x1;
    
        % 计算 
        temp=36*p1/hi*xi.^2.*(1-xi).^2+hi*p2*(1-xi).^4.*(2*xi+1).^2;                        a11=(5*temp(1)+8*temp(2)+5*temp(3))/9;
        temp=-6*p1*(xi-1).*(3*xi-1).*(1-xi).*xi+hi*hi*p2*xi.*(xi-1).^4.*(2*xi+1);           a12=(5*temp(1)+8*temp(2)+5*temp(3))/9;
        temp=-36*p1/hi*xi.^2.*(1-xi).^2.+hi*p2*xi.^2.*(1-xi).^2.*(3-2*xi).*(2*xi+1);        a13=(5*temp(1)+8*temp(2)+5*temp(3))/9;
        temp=-6*p1*(1-xi).*xi.*xi.*(3*xi-2)+hi*hi*p2*(xi-1).*xi.^2.*(1-xi).^2.*(2*xi+1);    a14=(5*temp(1)+8*temp(2)+5*temp(3))/9;
                                                                                            a21=a12;
        temp=hi*p1*(xi-1).^2.*(3*xi-1).^2+hi^3*p2*xi.^2.*(xi-1).^4;                         a22=(5*temp(1)+8*temp(2)+5*temp(3))/9;
        temp=6*p1*xi.*(1-xi).*(xi-1).*(3*xi-1)+hi*hi*p2*xi.^3.*(xi-1).^2.*(3-2*xi);         a23=(5*temp(1)+8*temp(2)+5*temp(3))/9;
        temp=p1*hi*(xi-1).*(3*xi-1).*xi.*(3*xi-2)+hi^3*p2*xi.^3.*(xi-1).^3;                 a24=(5*temp(1)+8*temp(2)+5*temp(3))/9;
                                                                                            a31=a13;
                                                                                            a32=a23;
        temp=36*p1/hi*xi.^2.*(1-xi).^2+hi*p2*xi.^4.*(3-2*xi).^2;                            a33=(5*temp(1)+8*temp(2)+5*temp(3))/9; 
        temp=6*p1*xi.*(3*xi-2).*xi.*(1-xi)+hi*hi*p2*(xi-1).*xi.^4.*(3-2*xi);                a34=(5*temp(1)+8*temp(2)+5*temp(3))/9;
                                                                                            a41=a14;
                                                                                            a42=a24;
                                                                                            a43=a34;
        temp=hi*p1*xi.^2.*(3*xi-2).^2+hi^3*p2*xi.^4.*(xi-1).^2;                             a44=(5*temp(1)+8*temp(2)+5*temp(3))/9;
    
        % 计算 b1, b2, b3, b4
        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*(1-xi).^2.*(1+2*xi) ; 
        b1=(5*temp(1)+8*temp(2)+5*temp(3))/9;

        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi))*hi.*xi.*(xi-1).^2;
        b2=(5*temp(1)+8*temp(2)+5*temp(3))/9; 
        
        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*xi.^2.*(3-2*xi);
        b3=(5*temp(1)+8*temp(2)+5*temp(3))/9;

        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*hi.*xi.^2.*(xi-1);
        b4=(5*temp(1)+8*temp(2)+5*temp(3))/9;
    
        % 组装刚度矩阵
        A(i1, i1)=A(i1, i1)+a11;  A(i1, i2)=A(i1, i2)+a12;  A(i1, i3)=A(i1, i3)+a13;  A(i1, i4)=A(i1, i4)+a14;
        A(i2, i1)=A(i2, i1)+a21;  A(i2, i2)=A(i2, i2)+a22;  A(i2, i3)=A(i2, i3)+a23;  A(i2, i4)=A(i2, i4)+a24;
        A(i3, i1)=A(i3, i1)+a31;  A(i3, i2)=A(i3, i2)+a32;  A(i3, i3)=A(i3, i3)+a33;  A(i3, i4)=A(i3, i4)+a34;
        A(i4, i1)=A(i4, i1)+a41;  A(i4, i2)=A(i4, i2)+a42;  A(i4, i3)=A(i4, i3)+a43;  A(i4, i4)=A(i4, i4)+a44;
        
        y(i1, 1)=y(i1, 1)+b1;     y(i2, 1)=y(i2, 1)+b2;     y(i3, 1)=y(i3, 1)+b3;     y(i4, 1)=y(i4, 1)+b4;
    end
    
    % 处理本质边界条件
    A(1, : )=0; A(1,1)=1;
 
    A(2*n+2, : )=0;A(2*n+2,2*n+2)=1;
    y(1)=0; y(2*n+2)=0;


    %% 求解系数矩阵
    u=zeros(2*n+2, 1);
    u=A\y;

    CondN(ni)=cond(A);

    %% 计算 L2 和 H1 误差
    % 利用复化Guass公式分区间段进行数值积分
    % 每个区间段内 u=u_{i1} N0(\xi)+ u_{ic} Nc(\xi)  + u_{i2} N1(\xi)
    % 每个区间段内 u'=u_{i1} N0'(\xi)+ u_{ic} Nc'(\xi)  + u_{i2} N1'(\xi)
    for i=1:1:n
        i1=I(1, i); i2=I(2, i); i3=I(3,i); i4=I(4,i);
        x1=p(i1);   x2=p(i3);   
        hi=x2-x1;

        % u 和 真解的值
        tempu=u(i1)*(1-xi).^2.*(2*xi+1)+ u(i2)*hi*xi.*(xi-1).^2+ u(i3)*xi.^2.*(3-2*xi)+ u(i4)*hi*xi.^2.*(xi-1);
        tempv=sin(pi/2*(x1+hi*xi));

        % u 和 真解的导数值
        tempdu=-6*u(i1)*xi.*(1-xi)/hi+u(i2)*(xi-1).*(3*xi-1)+u(i3)*6*xi.*(1-xi)/hi+u(i4)*xi.*(3*xi-2);
        tempdv=pi/2*cos(pi/2*(x1+hi*xi));

        % 计算 L2 误差
        temp=(tempu-tempv).^2*hi;
        err2(ni)=err2(ni)+(5*temp(1)+8*temp(2)+5*temp(3))/9;

        % 计算 H1 半范数
        temp=(tempdu-tempdv).^2*hi;
        errf(ni)=errf(ni)+(5*temp(1)+8*temp(2)+5*temp(3))/9;
    end

    % 计算 L2 和 H1 范数
    errH1(ni)=sqrt(errf(ni)+err2(ni));
    errL2(ni)=sqrt(err2(ni));
end

%% 画数值解

pp(1,[1:1:n+1])=a+[0:1:n]/n;
v=sin(pi/2*pp);
u=u([1:2:2*n+1],1);
figure
plot(pp,u,'-.r','DisplayName','u(x)');
title('数值解和精确解图像');
hold on
plot(pp,v,'--k','DisplayName','u*(x)');
legend('show');
hold off


%% 画误差图
% L^2 误差
N=10*[1:1:MK];
h=1./N;
figure
loglog(h, errL2, '--r', 'DisplayName', 'L^2 error');
hold on
loglog(h, h.^4, '-.k', 'DisplayName', 'O(h^4)')
legend('Show');
hold off

% H1 误差
figure
N=10*[1:1:MK];
h=1./N;
loglog(h, errH1, '--r', 'DisplayName', 'H1 error');
hold on
loglog(h, h.^3, '-.k', 'DisplayName', 'O(h^3)')
legend('Show');
hold off

%% Condition number
figure
loglog(h, CondN, '--r', 'DisplayName', 'Condition Number');
hold on
loglog(h, h.^-2, '-.k', 'DisplayName', 'y=1/h^{-2}');
legend('Show');
hold off


% 数值阶计算
alpha1=zeros(MK-1, 1);
alpha2=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errL2(i+1))-log(errL2(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(errH1(i+1))-log(errH1(i)))/(log(h(i+1))-log(h(i)));    
end

