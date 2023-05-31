%% 一维P3有限元直接生成刚度矩阵
%% 此程序为一维三次有限元方法， 基函数为Hermite 型基函数，节点为网格节点
%% 数值积分公式采用复化 Gauss 型求积公式

%% 测试数据为 a(u, v)=\int_{p u' v'+q uv}=(f, v) 于左端点满足本质边界条件
% p=1; q=pi^2/4
a=0; c=1;

MK=30;
errL2=zeros(MK, 1); errH=zeros(MK, 1); errH1=zeros(MK, 1);

for ni=1:1:MK
    %% Step 2: 生成网格信息
    N=10+10*ni;   % 区间剖分分数
    p=[0:1:N]*(c-a)/N;  % 以均匀网格为例

    % 声明总刚度矩阵
    A=zeros(2*N+2, 2*N+2);
    b=zeros(2*N+2, 1);
    
    % 收集网格信息
    G=zeros(6, N);
    for i=1:1:N
        G(1, i)=2*i-1;  G(2, i)=2*i;  G(3, i)=2*i+1; G(4, i)=2*i+2; 
        G(5, i)=p(i);   G(6, i)=p(i+1)-p(i);
    end
    
    %% Step 4: 形成刚度矩阵
    % 确定数值积分公式
    % 高斯型求积公式
    NM=4;
    xi = GaussianQPoints(0, 1, NM);

    % 计算基底向量 N(\xi)
    VecN=zeros(4, NM);    
    VecN(1, :)=(1-xi).^2.*(2*xi+1);
    VecN(2, :)=xi.*(1-xi).^2;
    VecN(3, :)=xi.^2.*(3-2*xi);
    VecN(4, :)=xi.^2.*(xi-1);

    % 计算基底向量 N'(\xi)
    VecM=zeros(4, NM);
    VecM(1, :)=-6*xi.*(1-xi);
    VecM(2, :)=(1-xi).*(1-3*xi);
    VecM(3, :)=6*xi.*(1-xi);
    VecM(4, :)=xi.*(3*xi-2); 

    % 此处采用不区分本质边界条件的方式
    for i=1:1:N
        % 计算单元刚度矩阵
        K=zeros(4, 4);  L=zeros(4, 1);
        x0=G(5, i); h=G(6, i);
    
        % 计算系数向量
        VecP=ones(1, NM);
        VecQ=pi^2/4*ones(1, NM);
        VecF=pi^2/2*sin(pi/2*(x0+h*xi));
    
        % 形成单元刚度矩阵
        for j1=1:1:4
            for j2=1:1:4
                temp=VecP.*VecM(j1, :).*VecM(j2,:)/h+VecQ.*VecN(j1, :).*VecN(j2, :)*h;   
                 K(j1, j2)=GaussianQuadrature(0, 1, temp, NM);
            end
            temp=VecF.*VecN(j1, :)*h;
            L(j1)=GaussianQuadrature(0, 1, temp, NM);
        end
    
        % 组装刚度矩阵
        for j1=1:1:4
            for j2=1:1:4
                A(G(j1, i), G(j2, i))=A(G(j1, i), G(j2, i))+K(j1, j2);
            end
            b(G(j1, i), 1)=b(G(j1, i), 1)+L(j1);
        end
    end
    
    % 处理本质边界条件
    A=A(2:1:2*N+2, 2:1:2*N+2);
    b=b(2:1:2*N+2, 1);
    
    %% Step 5: 求解系数矩阵
    u=A\b;
    u=[0; u];
    
    %% L2 误差和H1误差
    % u*=sin(pi/2*x);
    % u'*=pi/2cos(pi/2 x);
    for i=1:1:N
        x0=G(5, i); h=G(6, i);
    
        % 计算 L2 和 H1半范数在区间I(i) 上的值
        VecU=sin(pi/2*(x0+h*xi));
        VecDu=pi/2*cos(pi/2*(x0+h*xi));
    
        temp1=u(G(1, i))*VecN(1, :)+u(G(2, i))*VecN(2, :)+u(G(3, i))*VecN(3, :)+u(G(4, i))*VecN(4, :)-VecU;
        temp2=u(G(1, i))*VecM(1, :)/h+u(G(2, i))*VecM(2, :)/h+u(G(3, i))*VecM(3, :)/h+u(G(4, i))*VecM(4, :)/h-VecDu;

        errL2(ni)=errL2(ni)+h* GaussianQuadrature(0, 1, temp1.^2, NM);
        errH(ni)=errH(ni)+h*GaussianQuadrature(0, 1, temp2.^2, NM);

    end
    
    errH1(ni)=sqrt(errL2(ni)+errH(ni));
    errL2(ni)=sqrt(errL2(ni));
end

%% 误差图
%% 画误差图
% L^2 误差
N=10+10*[1:1:MK];
h=1./N;
loglog(h, errL2, '--r', 'DisplayName', 'L^2 error');
hold on
loglog(h, h.^4, '-.k', 'DisplayName', 'O(h^4)')
legend
hold off

% H1 误差
figure
N=10+10*[1:1:MK];
h=1./N;
loglog(h, errH1, '--r', 'DisplayName', 'H1 error');
hold on
loglog(h, h.^3, '-.k', 'DisplayName', 'O(h^3)')
legend
hold off

% 数值阶计算
alpha1=zeros(MK-1, 1);
alpha2=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errL2(i+1))-log(errL2(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(errH1(i+1))-log(errH1(i)))/(log(h(i+1))-log(h(i)));    
end

