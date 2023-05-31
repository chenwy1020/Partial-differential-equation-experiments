%% 二维矩形网格剖分乘积型Lagrange 双线性元
%% 数值积分公式采用乘积型 Gauss 型求积公式
%% 基于 sparse 形式

%% 测试问题 u*=xy+sin(pi x) sin(pi y)
% \Delta u+u=xy in [0, 1]*[0, 1]
% u(0, y)=u(x, 0)=0
% \partial_x u=y-\pi sin \pi y \in {x=1}*[0, 1]
% \partial_y u=x-\pi sin \pi x \in [0, 1]*{y=1}

%% 区域 \Omega=[a1, a2]*[a3, a4]
a1=0; a2=1; a3=0; a4=1;

NK=30;
L2err=zeros(NK, 1); H1err=zeros(NK, 1); H0err=zeros(NK, 1); meshsize=zeros(NK, 1);
CondN=zeros(NK, 1);
cpuN=zeros(NK, 1);

for ni=1:1:NK
    %% Step 2: 生成网格信息: p, I, E, G
    N=5*ni; M=5*ni;
    
    tic

    Rows=zeros(16*N*M+1, 1);
    Cols=zeros(16*N*M+1, 1);
    Vals=zeros(16*N*M+1, 1);

    % 声明总刚度矩阵
    b=zeros((N+1)*(M+1), 1);
    
    p=zeros(2, (N+1)*(M+1));
    I=zeros(4, N*M);
    E=zeros(3, 2*(N+M));
    
    % 依旧以均匀网格为例
    % 节点坐标
    hx=(a2-a1)/M;  hy=(a4-a3)/N;
    x=a1+[0:1:M]*hx;
    y=a3+[0:1:N]*hy;
    for i=1:1:N+1
        for j=1:1:M+1
            k=(i-1)*(M+1)+j;
            p(1, k)=x(j);   p(2, k)=y(i);
        end
    end
    
    % 矩形信息
    for i=1:1:N
        for j=1:1:M
            k=(i-1)*M+j;
            I(1, k)=(i-1)*(M+1)+j;  I(2, k)=(i-1)*(M+1)+j+1;
            I(3, k)=i*(M+1)+j;       I(4, k)=i*(M+1)+j+1;
        end
    end
    
    % 边界信息
    for i=1:1:M
        E(1, i)=i; E(2, i)=i+1; E(3, i)=1;
        E(1, i+M+N)=N*(M+1)+i;  E(2, i+M+N)=N*(M+1)+i+1;    E(3, i+M+N)=3;
    end
    
    for i=1:1:N
        E(1, i+M)=(M+1)*i; E(2, i+M)=(M+1)*(i+1); E(3, i+M)=2;
        E(1, i+2*M+N)=(M+1)*(i-1)+1;  E(2, i+2*M+N)=i*(M+1)+1;    E(3, i+2*M+N)=4;
    end
    
    % 广义坐标: 双线性元的广义坐标为网格坐标
    G=zeros(4, N*M);
    G=I;
    
    %% Step 3:  Lagrange 型双线性元
    % 确定数值积分公式
    % 对于二维求解节点写作 f(s, t)
    
    % 高斯型求积公式：此处仅给出 x, y 方向皆为 4 点的情形
    NMX=4; xi = GaussianQPoints(0, 1, NMX);
    NMY=4; eta=GaussianQPoints(0, 1, NMY);
    
    % 由于是乘积型基函数，生成一维标准基函数
    NX=[1-xi; xi];  DNX=[-ones(1, length(xi)); ones(1, length(xi))];
    NY=[1-eta; eta]; DNY=[-ones(1, length(eta)); ones(1, length(eta))];
    
    % 乘积型基底
    VecN=cell(1, 4); VecNx=cell(1, 4); VecNy=cell(1,4);
    for j1=1:1:2
        for j2=1:1:2
            % 节点基函数
            temp=NX(j2, :)'*NY(j1, :);  VecN{(j1-1)*2+j2}=temp;
    
            % 节点基函数\partial x
            temp=DNX(j2, :)'*NY(j1, :); VecNx{(j1-1)*2+j2}=temp;
    
            % 节点基函数\partial y
            temp=NX(j2, :)'*DNY(j1, :); VecNy{(j1-1)*2+j2}=temp;
        end
    end

    %% Step 4: 形成刚度矩阵
    % 此处采用不区分本质边界条件的方式
    % a(u, v)=\int_{\Omega} -q1ux vx-q2uy vy+puv...
    %   =int_{Omega} f1v-\int_{Gamma2}f2 vdy-\int_{Gamma3} f3 vdx
    for i=1:1:N*M
        % Step 4.2 计算单元刚度矩阵
        K=zeros(4,4);   L=zeros(4, 1);
        x0=p(1, G(1, i));   y0=p(2, G(1,i));
        hx=p(1, G(2, i))-x0;
        hy=p(2, G(3, i))-y0;
    
        % 定义参数矩阵 p, q1, q2
        VecP=2*pi^2*ones(NMX, NMY);
        VecQ1=-ones(NMX, NMY);
        VecQ2=-ones(NMX, NMY);
    
        % 定义右端函数 F
        VecF1=2*pi^2*(x0+hx*xi)'*(y0+hy*eta);
    
        % 形成单元刚度矩阵
        for j1=1:1:4
            for j2=1:1:4
                temp=VecP.*VecN{1, j1}.*VecN{1, j2}*hx*hy...
                    +VecQ1.*VecNx{1, j1}.*VecNx{1, j2}*hy/hx...
                    +VecQ2.*VecNy{1, j1}.*VecNy{1, j2}*hx/hy;
                K(j1, j2)=GaussQuadratureAB(0, 1, 0, 1, temp);
            end
            temp=VecF1.*VecN{1, j1}*hx*hy;
            L(j1)=GaussQuadratureAB(0, 1, 0, 1, temp);
        end
    
       %  Step 4.3 组装刚度矩阵
       for j1=1:1:4
           for j2=1:1:4
               Rows((j1-1)*4+j2+16*(i-1))=G(j1, i);
               Cols((j1-1)*4+j2+16*(i-1))=G(j2, i);
               Vals((j1-1)*4+j2+16*(i-1))=K(j1, j2);
           end
           b(G(j1, i), 1)=b(G(j1, i), 1)+L(j1);
       end
    end
    
    % Step 4.3 自然边界条件的处理
    % 自然边界条件：本问题中仅更新右端向量
    % Gamma2 边界条件的处理
    for i=1:1:N 
        x0=p(1, E(1, i+M));    y0=p(2, E(1, i+M));
        hy=p(2, E(2, i+M))-y0;
        VecF=-((y0+eta*hy)-pi*sin(pi*(y0+eta*hy)));
        L=zeros(2, 1);
    
        for j1=1:1:2
            temp=VecF.*NY(j1, :)*hy;
            L(j1)=GaussianQuadrature(0, 1, temp, NMY);
        end
    
        % 组装
        for j1=1:1:2
            b(E(j1, i+M), 1)= b(E(j1, i+M), 1)+L(j1);
        end
    end
    
    % Gamma3 边界条件的处理
    for i=1:1:M
        x0=p(1, E(1, i+M+N));   hx=p(1, E(2, i+M+N))-x0;
        VecF=-( (x0+hx*xi)-pi*sin(pi*(x0+hx*xi)));
        L=zeros(2, 1);
    
        for j1=1:1:2
            temp=VecF.*NX(j1, :)*hx;
            L(j1)=GaussianQuadrature(0, 1, temp, NMX);
        end
    
        % 组装
        for j1=1:1:2
            b(E(j1, i+M+N), 1)= b(E(j1, i+M+N), 1)+L(j1);
        end
    end
    
    Rows(16*N*M+1)=(N+1)*(M+1);
    Cols(16*N*M+1)=(N+1)*(M+1);
    Vals(16*N*M+1)=0;
    A=sparse(Rows, Cols, Vals);
    % 本质边界条件的处理，处理方式1：单位化相应行列
    temp=[E(1, 1:1:M), E(2, M), E(1, 1+2*M+N:1:2*M+2*N), E(2, 2*M+2*N)];
    A(temp, :)=0;  
    for i=1:1:length(temp)
        A(temp(i), temp(i))=1;
        b(temp(i))=0;
    end

    %% 条件数
    CondN(ni)=condest(A);
    
    %% Step 5: 系数矩阵求解
    u=A\b;


%     % 本质边界条件的处理，处理方式2：采用大数赋零
%     temp=[E(1, 1:1:M), E(2, M), E(1, 1+2*M+N:1:2*M+2*N), E(2, 2*M+2*N)];
%     for i=1:1:length(temp)
%         A(temp(i), temp(i))=10.^30;
%     end

%     %% Step 5: 系数矩阵求解
%     u=A\b;


    % 本质边界条件的处理，处理方式3：消去对应的行和列
%     temp=[E(1, 1:1:M), E(2, M), E(1, 1+2*M+N:1:2*M+2*N), E(2, 2*M+2*N)];
%     tempN=length(temp);
%     temp=sort(temp);  
%     % 消去在边界上的重复元
%     add1=1; temp1=temp(add1); 
%     for i=2:1:tempN
%         if(temp(i)==temp1)
%             temp(i)=0;
%         else
%             temp1=temp(i);
%             add1=add1+1;
%         end
%     end
%     temp=sort(temp);
%     temp=temp(tempN-add1+1:1:tempN);
% 
%     % 记录消去边界上的元之后的排序
%     Evec=zeros(1, add1); add2=1; add=0;
%     for i=1:1:(N+1)*(M+1)
%         if(add2<=add1)
%             if(i<temp(add2))
%                 add=add+1;
%                 Evec(1, add)=i;
%             else
%                 add2=add2+1;
%             end     
%         else
%             add=add+1;
%             Evec(1, add)=i;
%         end
%     end
% 
%     % 消去本质边界条件
%     A(temp, :)=[];  A(:, temp)=[]; b(temp)=[];
% 
%     %% Step 5: 系数矩阵求解
%     v=A\b;
% 
%     % 补回相应边界上的元
%     u=zeros((N+1)*(M+1), (N+1)*(M+1));
%     u(Evec)=v;  
    
    %% Step 6: 误差计算
    for i=1:1:N*M
        x0=p(1, G(1, i));   y0=p(2, G(1,i));
        hx=p(1, G(2, i))-x0;
        hy=p(2, G(3, i))-y0;
    
        % 真解及其偏导数
        tempx=x0+hx*xi;     tempy=y0+hy*eta;
        utrue=tempx'*tempy+sin(pi*tempx)'*sin(pi*tempy);
        dxu=ones(NMX, 1)*tempy+pi*cos(pi*tempx)'*sin(pi*tempy);
        dyu=tempx'*ones(1, NMY)+pi*sin(pi*tempx)'*cos(pi*tempy);
    
        tempu=zeros(NMX, NMY);
        tempdxu=zeros(NMX, NMY);
        tempdyu=zeros(NMX, NMY);
        for j=1:1:4
                tempu=tempu+u(G(j, i))*VecN{1, j};
                tempdxu=tempdxu+u(G(j, i))*VecNx{1, j}/hx;
                tempdyu=tempdyu+u(G(j, i))*VecNy{1, j}/hy;
        end
        tempu=tempu-utrue;
        tempdxu=tempdxu-dxu;
        tempdyu=tempdyu-dyu;
    
        L2err(ni)=L2err(ni)+hx*hy*GaussQuadratureAB(0, 1, 0, 1, tempu.^2);
        H0err(ni)=H0err(ni)+hx*hy*GaussQuadratureAB(0, 1, 0, 1, tempdxu.^2)...
            +hx*hy*GaussQuadratureAB(0, 1, 0, 1, tempdyu.^2);
    end

H1err(ni)=L2err(ni)+H0err(ni);
L2err(ni)=sqrt(L2err(ni));
H1err(ni)=sqrt(H1err(ni));

% 网格最大边长
meshsize(ni)=sqrt(hx^2+hy^2);

toc
cpuN(ni)=toc;

end

%% 作图
loglog(meshsize, L2err, '-*', 'DisplayName', 'L2err');
hold on
loglog(meshsize, meshsize.^2, '--k', 'DisplayName', 'O(h^2)')
legend('Show')
hold off

figure
loglog(meshsize, H1err, '-*', 'DisplayName', 'H1err');
hold on
loglog(meshsize, meshsize, '--k', 'DisplayName', 'O(h)')
legend('Show')
hold off

figure
loglog(meshsize, CondN, '-*', 'DisplayName', 'Condition number');
hold on
loglog(meshsize, meshsize.^-2, '--k', 'DisplayName', 'O(1/h^2)')
legend('Show')
hold off

figure
loglog(meshsize, cpuN, '-*', 'DisplayName', 'Complexity');
hold on
loglog(meshsize, meshsize.^-2, '--k', 'DisplayName', 'O(1/h^2)')
legend('Show')
hold off

%% 数值收敛阶
alpha1=zeros(NK-1, 1);
alpha2=zeros(NK-1, 1);
for i=1:1:NK-1
    alpha1(i)=(log(L2err(i+1))-log(L2err(i)))/(log(meshsize(i+1))-log(meshsize(i)));
    alpha2(i)=(log(H1err(i+1))-log(H1err(i)))/(log(meshsize(i+1))-log(meshsize(i)));    
end


















