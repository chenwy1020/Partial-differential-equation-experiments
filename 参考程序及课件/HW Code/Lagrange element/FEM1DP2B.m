%% 一维P2有限元去掉中间节点
%% 此程序为一维二次有限元方法， 基函数为Lagrange型基函数，节点为网格节点和中点
%% 数值积分公式采用复化Simpson公式和Gauss型求积公式以及复化梯形公式

%% 测试数据为 a(u, v)=\int_{p u' v'+q uv}=(f, v) 于左端点满足本质边界条件
% p=1; q=pi^2/4
a=0; c=1;

MK=50;
errL2=zeros(MK, 1); errH=zeros(MK, 1); errH1=zeros(MK, 1);
cpus=zeros(MK, 1);

for ni=1:1:MK
    tic;
    %% Step 2: 生成网格信息
    N=10+2*ni;   % 区间剖分分数
    p=[0:1:N]*(c-a)/N;  % 以均匀网格为例
    
    % 声明总刚度矩阵：仅含有节点
    A=zeros(N+1, N+1);
    b=zeros(N+1, 1);
    
    % 收集网格信息
    G=zeros(4, N);
    for i=1:1:N
        G(1, i)=i;  G(2, i)=i+1;   
        G(3, i)=p(i);   G(4, i)=p(i+1)-p(i);
    end
    
    % 用于存储中点相关的权重
    HalfP=zeros(3, N);
    
    %% Step 4: 形成刚度矩阵
    % 确定数值积分公式
    % 复化Simpson公式
    M=4; hm=1/(2*M); xi=[0:1:2*M]*hm; NM=2*M+1;

%     %高斯型求积公式
%     NM=4;
%     xi = GaussianQPoints(0, 1, NM);

%     % 复化梯形求积公式
%     NM=10; xi=[0:1:NM-1]*1/(NM-1);

    % 计算基底向量 N(\xi)
    VecN=zeros(3, NM);    
    VecN(1, :)=(2*xi-1).*(xi-1);
    VecN(2, :)=(2*xi-1).*xi;
    VecN(3, :)=4*xi.*(1-xi);

    % 计算基底向量 N'(\xi)
    VecM=zeros(3, NM);
    VecM(1, :)=(4*xi-3);
    VecM(2, :)=(4*xi-1);
    VecM(3, :)=(4-8*xi);

    % 此处采用不区分本质边界条件的方式
    for i=1:1:N
        % 计算单元刚度矩阵
        K=zeros(3, 3);  L=zeros(3, 1);
        x0=G(3, i); h=G(4, i);
    
        % 计算系数向量
        VecP=ones(1, NM);
        VecQ=pi^2/4*ones(1, NM);
        VecF=pi^2/2*sin(pi/2*(x0+h*xi));
    
        % 行程单元刚度矩阵
        for j1=1:1:3
            for j2=1:1:3
                temp=VecP.*VecM(j1, :).*VecM(j2,:)/h+VecQ.*VecN(j1, :).*VecN(j2, :)*h;
                K(j1, j2)=QuadratureRules(0, 1, temp, 1);
            end
            temp=VecF.*VecN(j1, :)*h;
            L(j1)=QuadratureRules(0, 1, temp, 1);
        end
    
        %消去中间节点
        HalfP(1, i)=-K(1, 3)/K(3, 3);
        HalfP(2, i)=-K(2, 3)/K(3, 3);
        HalfP(3, i)=L(3)/K(3, 3);
    
        for j1=1:1:2
            for j2=1:1:2
                K(j1, j2)=K(j1, j2)-K(3, j2)/K(3, 3)*K(j1, 3);
            end
            L(j1)=L(j1)-K(j1, 3)/K(3,3)*L(3);
        end
    
        % 组装刚度矩阵   
        for j1=1:1:2
                for j2=1:1:2
                    A(G(j1, i), G(j2, i))=A(G(j1, i), G(j2, i))+K(j1, j2);
                end
                b(G(j1, i), 1)=b(G(j1, i), 1)+L(j1);
        end
    end
    
    % 处理本质边界条件
    A=A(2:1:N+1, 2:1:N+1);
    b=b(2:1:N+1, 1);
    
    %% Step 5: 求解系数矩阵
    u=A\b;
    u=[0; u];
    
    % 获得中点值
    halfu=HalfP(1, :).*u(1:1:N, 1)'+HalfP(2, :).*u(2:1:N+1, 1)'+HalfP(3, :);

  %% L2 误差和H1误差
    % u*=sin(pi/2*x);
    % u'*=pi/2cos(pi/2 x);
    for i=1:1:N
        x0=G(3, i); h=G(4, i);
    
        % 计算 L2 和 H1半范数在区间I(i) 上的值
        VecU=sin(pi/2*(x0+h*xi));
        VecDu=pi/2*cos(pi/2*(x0+h*xi));
    
        temp1=u(G(1, i))*VecN(1, :)+u(G(2, i))*VecN(2, :)+halfu(i)*VecN(3, :)-VecU;
        temp2=u(G(1, i))*VecM(1, :)/h+u(G(2, i))*VecM(2, :)/h+halfu(i)*VecM(3, :)/h-VecDu;
    
        errL2(ni)=errL2(ni)+h*QuadratureRules(0, 1, temp1.^2, 1);
        errH(ni)=errH(ni)+h*QuadratureRules(0, 1, temp2.^2, 1);
    end
    
    errH1(ni)=sqrt(errL2(ni)+errH(ni));
    errL2(ni)=sqrt(errL2(ni));

    cpus(ni)=toc;
end

%% 误差图
%% 画误差图
% L^2 误差
N=10+2*[1:1:MK];
h=1./N;
loglog(h, errL2, '--r', 'DisplayName', 'L^2 error');
hold on
loglog(h, h.^3, '-.k', 'DisplayName', 'O(h^3)')
legend
hold off

% H1 误差
figure
N=10+2*[1:1:MK];
h=1./N;
loglog(h, errH1, '--r', 'DisplayName', 'H1 error');
hold on
loglog(h, h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend
hold off

% 数值阶计算
alpha1=zeros(MK-1, 1);
alpha2=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errL2(i+1))-log(errL2(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(errH1(i+1))-log(errH1(i)))/(log(h(i+1))-log(h(i)));    
end








