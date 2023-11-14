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
    A=zeros(2*n+1, 2*n+1); y=zeros(2*n+1, 1);
    p=zeros(1, 2*n+1); I=zeros(3, n);
    
    %% 网格信息
    % 这里考虑均匀剖分 
    h=(b-a)/(2*n);
    p(1, : )=a+[0:1:2*n]*h; 
    I(1, : )=(1:2:2*n-1);  I(2, : )=(3:2:2*n+1); I(3, : )=(2:2:2*n);

    %% 给出变分问题 a(u, v)
    % 对于此问题 p(x)=1, q(x)=pi^2/4
    p1=1; p2=pi^2/4;

    %% 定义标准区域上的数值求积公式，定义求积节点
    % 这里利用复化Guass公式
    xi(1)=(a+b)/2-(b-a)/2*sqrt(3/5);
    xi(2)=(a+b)/2;
    xi(3)=(a+b)/2+(b-a)/2*sqrt(3/5);
    
    %% 计算单元刚度矩阵
    for i=1:1:n
        i1=I(1, i); i2=I(2, i); ic=I(3,i);
        x1=p(i1);  x2=p(i2);  xc=p(ic);
        hi=x2-x1;
    
        % 计算 a1, a2, a3, a4
        temp=p1/hi*(4*xi-3).^2+hi*p2*(2*xi-1).^2.*(xi-1).^2;                        a11=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
        temp=4*p1/hi*(4*xi-3).*(1-2*xi)+4*hi*p2*(2*xi-1).*(xi-1).*xi.*(1-xi);       a12=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
        temp=p1/hi*(4*xi-1).*(4*xi-3)+hi*p2*(2*xi-1).^2.*xi.*(xi-1);                a13=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
                                                                                    a21=a12;
        temp=16*p1/hi*(1-2*xi).^2+16*hi*p2*(1-xi).^2.*xi.^2;                        a22=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
        temp=-4*p1/hi*(1-4*xi).*(1-2*xi)+4*hi*p2*(2*xi-1).*xi.^2.*(1-xi);           a23=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
                                                                                    a31=a13;
                                                                                    a32=a23;
        temp=p1/hi*(4*xi-1).^2+hi*p2*(2*xi-1).^2.*xi.^2;                            a33=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0; %可以和a11相加
    
        % 计算 b1, b2
        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*(1-2*xi).*(-xi);%左侧积分 实际上应当是右侧，不过首项系数为零，也就罢了
        b1=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;

        temp=4*hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*xi.*(1-xi);
        bc=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0; 

        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*(2*xi-1).*(xi-1);%右侧积分
        b2=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
    
        % 组装刚度矩阵
        A(i1, i1)=A(i1, i1)+a11;  A(i1, ic)=A(i1, ic)+a12;  A(i1, i2)=A(i1, i2)+a13;
        A(ic, i1)=A(i2, i2)+a21;  A(ic, ic)=A(ic, ic)+a22;  A(ic, i2)=A(ic, i2)+a23;
        A(i2, i1)=A(i2, i1)+a31;  A(i2, ic)=A(i2, ic)+a32;  A(i2, i2)=A(i2, i2)+a33;

        y(i1, 1)=y(i1, 1)+b1;     y(ic, 1)=y(ic, 1)+bc;     y(i2, 1)=y(i2, 1)+b2;
    end
    
    % 处理本质边界条件
    A=A(2:1:2*n+1, 2:1:2*n+1);
    y=y(2:1:2*n+1);

    %% 求解系数矩阵
    u=zeros(2*n, 1);
    u=A\y;
    u=[0; u];
    
    CondN(ni)=cond(A);

    %% 计算 L2 和 H1 误差
    % 利用复化Guass公式分区间段进行数值积分
    % 每个区间段内 u=u_{i1} N0(\xi)+ u_{ic} Nc(\xi)  + u_{i2} N1(\xi)
    % 每个区间段内 u'=u_{i1} N0'(\xi)+ u_{ic} Nc'(\xi)  + u_{i2} N1'(\xi)
    for i=1:1:n
        i1=I(1, i); i2=I(2, i); ic=I(3,i);
        x1=p(i1);   x2=p(i2);   xc=p(ic);
        hi=x2-x1;

        % u 和 真解的值
        tempu=u(i1)*(2*xi-1).*(xi-1)+u(ic)*4*xi.*(1-xi)+u(i2)*(1-2*xi).*(-xi);
        tempv=sin(pi/2*(x1+hi*xi));

        % u 和 真解的导数值
        tempdu=u(i1)*(4*xi-3)/hi+u(ic)*4*(1-2*xi)/hi+u(i2)*(4*xi-1)/hi;
        tempdv=pi/2*cos(pi/2*(x1+hi*xi));

        % 计算 L2 误差
        temp=(tempu-tempv).^2*hi/2.0;
        err2(ni)=err2(ni)+(5*temp(1)+8*temp(2)+5*temp(3))/9;

        % 计算 H1 半范数
        temp=(tempdu-tempdv).^2*hi/2.0;
        errf(ni)=errf(ni)+(5*temp(1)+8*temp(2)+5*temp(3))/9;
    end

    % 计算 L2 和 H1 范数
    errH1(ni)=sqrt(errf(ni)+err2(ni));
    errL2(ni)=sqrt(err2(ni));
end

%% 画数值解

v=sin(pi/2*p);

figure
plot(p,u,'-.r','DisplayName','u(x)');
title('数值解和精确解图像');
hold on
plot(p,v,'--k','DisplayName','u*(x)');
legend('show');
hold off


%% 画误差图
% L^2 误差
N=10*[1:1:MK];
h=1./N;
figure
loglog(h, errL2, '--r', 'DisplayName', 'L^2 error');
hold on
loglog(h, h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

% H1 误差
figure
N=10*[1:1:MK];
h=1./N;
loglog(h, errH1, '--r', 'DisplayName', 'H1 error');
hold on
loglog(h, h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

%% Condition number
figure
loglog(h, CondN, '--r', 'DisplayName', 'Condition Number');
hold on
loglog(h, h.^-2, '-.k', 'DisplayName', 'y=h^{-2}');
legend('Show');
hold off


% 数值阶计算
alpha1=zeros(MK-1, 1);
alpha2=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errL2(i+1))-log(errL2(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(errH1(i+1))-log(errH1(i)))/(log(h(i+1))-log(h(i)));    
end

