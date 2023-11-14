%此程序为Lagrange型一维二次有限元有限元方法的通用程序
clear
clc

% 进行循环的次数
MK=6;

err2=zeros(MK, 1);      errf=zeros(MK, 1);      
errH1=zeros(MK, 1);   errL2=zeros(MK,1);
CondN=zeros(MK, 1);
p=@(x)(1);
xmin=0; xmax=1;

for ni=1:1:MK
    n=10*ni;
    A=zeros(2*n+1, 2*n+1);   u=zeros(2*n+1, 1);     y=zeros(2*n+1, 1);     
    points=zeros(2*n+1, 1);  dualpoints=zeros(2*n+2, 1);
   
    %% 网格信息
    % 均匀剖分,均匀对偶剖分
    h=(xmax-xmin)/n;
    points(:,1)=xmin+(0:1:2*n)*h/2.0;   
    dualpoints(2:1:2*n+1,1)=xmin+h/4.0+(0:1:2*n-1)*h/2.0;
    dualpoints(1)=xmin;    dualpoints(2*n+2)=xmax;
    
    %% 计算系数矩阵和右端向量
    % i是遍历对偶区间的指标（设[x0,x1/4]为第2项）
    for i=2:1:2*n
        if mod(i,2)==0
            %% 半结点points(i) 对应的对偶区间
            % 左右整结点
            x1=points(i-1); x2=points(i+1);
            hi=x2-x1;
            % 左右对偶结点
            xh1=dualpoints(i); xh2=dualpoints(i+1);

            % 定义标准区域上的数值求积公式，定义求积节点，这里利用复化Guass公式
            a=(xh1-x1)/hi; b=(xh2-x1)/hi;

            % 计算 a1, a2, a3, b1
            a1=p(xh1)*( 4*a-3)/hi -p(xh2)*( 4*b-3)/hi +Gauss(a,b,x1,hi,1);
            a2=p(xh1)*(-8*a+4)/hi -p(xh2)*(-8*b+4)/hi +Gauss(a,b,x1,hi,2);
            a3=p(xh1)*( 4*a-1)/hi -p(xh2)*( 4*b-1)/hi +Gauss(a,b,x1,hi,3);
            b1=Gauss(a,b,x1,hi,0);

            % 组装
            A(i,i-1)=a1;
            A(i,i  )=a2;
            A(i,i+1)=a3;
            y(i)=b1;
        else
            %% 整结点points(i) 对应的对偶区间
            % 左右整结点
            x1=points(i-2); x2=points(i); x3=points(i+2);
            hi1=x2-x1; hi2=x3-x2;

            % 左右对偶结点
            xh1=dualpoints(i); xh2=dualpoints(i+1);

            % 计算 a1, a2, a3, a4, a5, b1
            a=(xh1-x1)/hi1; b=1;
            a1=p(xh1)*( 4*a-3)/hi1+Gauss(a,b,x1,hi1,1);
            a2=p(xh1)*(-8*a+4)/hi1+Gauss(a,b,x1,hi1,2); 
            a3=p(xh1)*( 4*a-1)/hi1+Gauss(a,b,x1,hi1,3);
            b1=Gauss(a,b,x1,hi1,0);

            a=0; b=(xh2-x2)/hi2;
            a3=-p(xh2)*( 4*b-3)/hi2+Gauss(a,b,x2,hi2,1)+a3;
            a4=-p(xh2)*(-8*b+4)/hi2+Gauss(a,b,x2,hi2,2);
            a5=-p(xh2)*( 4*b-1)/hi2+Gauss(a,b,x2,hi2,3);
            b1=Gauss(a,b,x2,hi2,0)+b1;

            % 组装
            A(i,i-2)=a1;
            A(i,i-1)=a2;
            A(i,i  )=a3;
            A(i,i+1)=a4;
            A(i,i+2)=a5;
            y(i)=b1;

        end
    end
    
    %% 初边值条件处理条件
    % i=1
    A(1,1)=1; y(1)=0;
    % i=2*n+1
    % 左右整结点
    i=2*n+1;
    x1=points(i-2); x2=points(i); 
    hi1=x2-x1; 

    % 左右对偶结点
    xh1=dualpoints(i);

    % 计算 a1, a2, a3, a4, a5, b1
    a=(xh1-x1)/hi1; b=1;
    a1=p(xh1)*( 4*a-3)/hi1+Gauss(a,b,x1,hi1,1);
    a2=p(xh1)*(-8*a+4)/hi1+Gauss(a,b,x1,hi1,2);
    a3=p(xh1)*( 4*a-1)/hi1+Gauss(a,b,x1,hi1,3);
    b1=Gauss(a,b,x1,hi1,0);

    % 组装
    A(2*n+1,2*n-1)=a1;
    A(2*n+1,2*n  )=a2;
    A(2*n+1,2*n+1)=a3;
    y(2*n+1)=b1;

    %A=sparse(A);

    %% 求解系数矩阵
    
    u=A\y;

    CondN(ni)=condest(A);
    
    %% 计算 L2 和 H1 误差
    % 利用复化Guass公式分区间段进行数值积分
    a=0; b=1;
    xi(1)=(a+b)/2.0-(b-a)/2.0*sqrt(3/5.0);
    xi(2)=(a+b)/2.0;
    xi(3)=(a+b)/2.0+(b-a)/2.0*sqrt(3/5.0);

    for i=1:1:n
        % u 和 真解的值
        x1=points(2*i-1); x2=points(2*i+1);
        hi=x2-x1;

        tempu=u(2*i-1)*(2*xi-1).*(xi-1)+u(2*i)*4*xi.*(1-xi)+u(2*i+1)*(2*xi-1).*xi;
        tempv=sin(pi/2*(x1+hi*xi));

        % u 和 真解的导数值
        tempdu=u(2*i-1)*(4*xi-3)/hi+u(2*i)*4*(1-2*xi)/hi+u(2*i+1)*(4*xi-1)/hi;
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

v=sin(pi/2*points);

figure
plot(points,u,'-.r','DisplayName','u(x)');
title('数值解和精确解图像');
hold on
plot(points,v,'--k','DisplayName','u*(x)');
legend('show');
hold off

%% 画误差图
% L^2 误差
N=10*(1:1:MK);
h=1./N;
figure
loglog(h, errL2, '--r', 'DisplayName', 'L^2 error');
hold on
loglog(h, h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

% H1 误差
figure
N=10*(1:1:MK);
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
% 数值阶计算
alpha1=zeros(MK-1, 1);
alpha2=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errL2(i+1))-log(errL2(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(errH1(i+1))-log(errH1(i)))/(log(h(i+1))-log(h(i)));    
end

