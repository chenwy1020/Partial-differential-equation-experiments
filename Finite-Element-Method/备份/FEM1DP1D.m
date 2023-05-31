%% 此程序为一维一次有限元方法的通用程序
clear

% 进行循环的次数
MK=100;

err2=zeros(MK, 1);      errf=zeros(MK, 1);      
errH1=zeros(MK, 1);   errL2=zeros(MK,1);
CondN=zeros(MK, 1);

for ni=1:1:MK
    n=10+2*ni;
    a=0; b=1;
    A=zeros(n+1, n+1); y=zeros(n+1, 1);
    p=zeros(1, n+1); I=zeros(2, n);
    
    %% 网格信息
    % 这里考虑均匀剖分 
    h=(b-a)/n;
    p=a+[0:1:n]*h;
    I(1, :)=[1:1:n];    I(2, : )=[2:1:n+1];
    
    %% 给出变分问题 a(u, v)
    % 对于此问题 p(x)=1, q(x)=pi^2/4
    p1=1; p2=pi^2/4;

    %% 定义标准区域上的数值求积公式，定义求积节点
    % 这里利用复化Simpson公式
    M=2; hm=1/(2*M); xi=[0:1:2*M]*hm;
    
    %% 计算单元刚度矩阵
    for i=1:1:n
        i1=I(1, i); i2=I(2, i);
        x1=p(i1);  x2=p(i2); hi=x2-x1;
    
        % 计算 a1, a2, a3, a4
        temp=-p1/hi+hi*p2*(1-xi).*xi;     a1=ComSimpson(0, 1, temp);
        temp=-p1/hi+hi*p2*(1-xi).*xi;     a2=ComSimpson(0, 1, temp);
        temp=p1/hi+hi*p2*xi.^2;   a3=ComSimpson(0, 1, temp);
        temp=p1/hi+hi*p2*(1-xi).^2;   a4=ComSimpson(0, 1, temp);
    
        % 计算 b1, b2
        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*(1-xi);
        b1=ComSimpson(0, 1, temp);
    
        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*xi;
        b2=ComSimpson(0, 1, temp);
    
        % 组装刚度矩阵
        A(i1, i2)=A(i1, i2)+a1;     A(i2, i1)=A(i2, i1)+a2;
        A(i2, i2)=A(i2, i2)+a3;     A(i1, i1)=A(i1, i1)+a4;
        y(i1, 1)=y(i1, 1)+b1;           y(i2, 1)=y(i2, 1)+b2;
    end
    
    % 处理本质边界条件：第一行赋值 [1, 0,...]
    A(1, 1)=10^10; 
    y(1)=0;

    %% 求解系数矩阵
    u=zeros(n+1, 1);
    u=A\y;

    CondN(ni)=cond(A);
    
    %% 计算 L2 和 H1 误差
    % 利用复化Simpson公式分区间段进行数值积分
    % 每个区间段内 u=u_{i1} N0(\xi)+u_{i2}N1(\xi)
    % 每个区间段内 u'=(u_{i2}-u_{i1})/h
    for i=1:1:n
        i1=I(1, i); i2=I(2, i);
        x1=p(i1);  x2=p(i2); hi=x2-x1;

        % u 和 真解的值
        tempu=u(i1)*(1-xi)+u(i2)*xi;
        tempv=sin(pi/2*(x1+hi*xi));

        % u 和 真解的导数值
        tempdu=(u(i2)-u(i1))/hi;
        tempdv=pi/2*cos(pi/2*(x1+hi*xi));

        % 计算 L2 误差
        temp=(tempu-tempv).^2*hi;
        err2(ni)=err2(ni)+ComSimpson(0, 1, temp);

        % 计算 H1 半范数
        temp=(tempdu-tempdv).^2*hi;
        errf(ni)=errf(ni)+ComSimpson(0, 1, temp);
    end

    % 计算 L2 和 H1 范数
    errH1(ni)=sqrt(errf(ni)+err2(ni));
    errL2(ni)=sqrt(err2(ni));
end

%% 画误差图
% L^2 误差
N=10+2*[1:1:MK];
h=1./N;
loglog(h, errL2, '--r', 'DisplayName', 'L^2 error');
hold on
loglog(h, h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

% H1 误差
figure
N=10+2*[1:1:MK];
h=1./N;
loglog(h, errH1, '--r', 'DisplayName', 'H1 error');
hold on
loglog(h, h*10, '-.k', 'DisplayName', 'O(h)')
legend('Show');
hold off

%% Condition number
figure
loglog(h, CondN, '--r', 'DisplayName', 'Condition Number');
hold on
loglog(h, h.^-2, '-.k', 'DisplayName', 'y=1/h^2');
legend('Show');
hold off


% 数值阶计算
alpha1=zeros(MK-1, 1);
alpha2=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errL2(i+1))-log(errL2(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(errH1(i+1))-log(errH1(i)))/(log(h(i+1))-log(h(i)));    
end

