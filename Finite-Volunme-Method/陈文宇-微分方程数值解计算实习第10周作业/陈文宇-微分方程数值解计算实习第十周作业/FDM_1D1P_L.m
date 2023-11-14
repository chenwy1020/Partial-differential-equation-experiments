 %仿照老师所给的程序代码，建立Lagrange型二次有限元方法
%此程序为Lagrange型一维二次有限元有限元方法的通用程序
clear
clc

% 进行循环的次数
MK=20;

errC=zeros(MK, 1);    err0=zeros(MK, 1);      err1=zeros(MK, 1); 
CondN=zeros(MK, 1);
p=@(x)(ones(length(x),1));
r=@(x)(zeros(length(x),1));
q=@(x)(ones(length(x),1)*pi^2/4);
f=@(x)(sin(pi*x/2)*pi^2/2);
vf=@(x)(sin(pi*x/2));


for ni=1:1:MK
    n=10*ni;
    a=0; b=1;
    A=zeros(n+1, n+1); y=zeros(n+1, 1);     u=zeros(n+1, 1);
    ph=zeros(n+1, 1);  pc=zeros(n+1, 1);
    I=zeros(n,3);
    
    %% 网格信息
    % 这里考虑均匀剖分 
    h=(b-a)/n;
    ph(:,1)=a+(0:1:n)*h;    pc(:,1)=a+(0:1:n)*h+h/2;
    I(:,1)=(1:1:n);         I(:,2)=(2:1:n+1);         I(:,3)=(3:1:n+2);

    %% 给出相关节点处的各类函数值
    
    VecP=p(pc);  VecR=r(ph);  VecQ=q(ph);  VecF=f(ph); Vecvf=vf(ph);
    
    %% 计算系数矩阵和右端向量
    for i=1:1:n-1
        i1=I(i,1); i2=I(i,2); i3=I(i,3);
        x1=ph(i1); x2=ph(i2); x3=ph(i3);
        hi1=x2-x1; hi2=x3-x2;
    
        % 计算 a1, a2, a3, b1
        %直接差分法 %有限体积法的中矩形公式
        a1=-VecP(i1)/hi1;
        a2=VecP(i2)/hi2+VecP(i1)/hi1+VecQ(i2)*(hi1+hi2)/2;
        a3=-VecP(i2)/hi2;
        b1=VecF(i2)*(hi1+hi2)/2;
        

        % 组装
        A(i2,i1)=a1;
        A(i2,i2)=a2;
        A(i2,i3)=a3;
        y(i2)=b1;

    end

    % 处理初边值条件
    
    %直接差分法
    %向后差分法
    % A(1,1)=1;  A(n+1,n)=-1/hi1;  A(n+1,n+1)=1/hi2;
    % y(1)=0;    y(n+1)=0;
    %中心差分
    % A(1,1)=1;  A(n+1,n)=-VecP(i2)/hi2;  A(n+1,n+1)=VecP(i2)/hi2+VecQ(i3)*hi2/2;
    % y(1)=0;    y(n+1)=VecF(i3)*hi2/2;

    %有限体积法的中矩形公式
    A(1,1)=1;  A(n+1,n)=-VecP(i2)/hi2;  A(n+1,n+1)=VecP(i2)/hi2+VecQ(i3)*hi2/2;
    y(1)=0;    y(n+1)=VecF(i3)*hi2/2;

    
    A=sparse(A);

    %% 求解系数矩阵
    
    u=A\y;

    CondN(ni)=condest(A);

    %% 计算误差
    temp=u(2:1:n)-Vecvf(2:1:n);
    errC(ni)=norm(temp,inf);
    
    temp=u(2:1:n)-Vecvf(2:1:n);
    err0(ni)=sqrt(norm(temp,2)^2 * h);
    
    temp=(u(2:1:n+1)-u(1:1:n))/h-(Vecvf(2:1:n+1)-Vecvf(1:1:n))/h;
    err1(ni)=sqrt(norm(temp,2)^2 * h)+err0(ni);

end

%% 画数值解

v=sin(pi/2*ph);

figure
plot(ph,u,'-.r','DisplayName','u(x)');
title('数值解和精确解图像');
hold on
plot(ph,v,'--k','DisplayName','u*(x)');
legend('show');
hold off


%% 画误差图
% norm-C 误差
N=10*(1:1:MK);
h=1./N;
figure
loglog(h, errC, '--r', 'DisplayName', 'errC');
hold on
loglog(h, h.^2/2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

% norm-0 误差
N=10*(1:1:MK);
h=1./N;
figure
loglog(h, err0, '--r', 'DisplayName', 'err0');
hold on
loglog(h, h.^2/2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

% norm-1 误差
figure
N=10*(1:1:MK);
h=1./N;
loglog(h, err1, '--r', 'DisplayName', 'err1');
hold on
loglog(h, h.^2/20, '-.k', 'DisplayName', 'O(h^{2})')
legend('Show');
hold off

%% Condition number
figure
loglog(h, CondN, '--r', 'DisplayName', 'Condition Number');
hold on
loglog(h, h.^-2/4, '-.k', 'DisplayName', 'y=h^{-2}');
legend('Show');
hold off


% 数值阶计算
alphaC=zeros(MK-1, 1);
alpha0=zeros(MK-1, 1);
alpha1=zeros(MK-1, 1);
for i=1:1:MK-1
    alphaC(i)=(log(errC(i+1))-log(errC(i)))/(log(h(i+1))-log(h(i)));
    alpha0(i)=(log(err0(i+1))-log(err0(i)))/(log(h(i+1))-log(h(i)));    
    alpha1(i)=(log(err1(i+1))-log(err1(i)))/(log(h(i+1))-log(h(i)));  
end

