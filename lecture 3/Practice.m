%% practice 2: 例1.4.1
%% 取基底 phi_i=sin(i pi x) 和 phi_i=(1-x)x^i 计算例 1.4.1

M=50;

err1=zeros(1,M);
err2=zeros(1,M);

condA=zeros(1, M);
condB=zeros(1, M);

x=[0:1:1000]/1000;
u=sin(x)/sin(1)-x;

%% 基底1
for i=1:1:M
    N=i;

    %% 初始化
    A=zeros(N, N);
    b=zeros(N, 1);

    %% 仅对角线元素不为零
    for k=1:1:N
        A(k, k)=((k*pi)^2-1)/2;
        b(k)=-cos(k*pi)/(k*pi);
    end

    %% 线性方程组求解
    S=A\b;
    
    %% 计算L2 误差
    % 计算 un 在求积节点处取值 
    C=zeros(1,1001);
    for k=1:1:N
        C=C+S(k)*sin(k*pi*x);
    end
    E=(C-u).^2;

    % 利用复化Simpson公式计算L^2 误差
    err1(i)=sqrt(ComSimpson(0, 1, E));
    
    if(i==10) 
       figure
       spy(A);
    end

    %% 计算条件数，cond(A, 2) 计算 A的从属 2 范数
    condA(i)=cond(A, 2);
end



%% 基底2
for i=1:1:M
    N=i;
    A=zeros(N, N);
    b=zeros(N, 1);

    %% 生成系数矩阵
    for k1=1:1:N
        for k2=1:1:N
            A(k1, k2)=k1*k2/(k1+k2-1)-(2*k1*k2+k1+k2)/(k1+k2)+(k1*k2+k1+k2)/(k1+k2+1)...
                +2/(k1+k2+2)-1/(k1+k2+3);
        end
        b(k1)=1/(k1+2)-1/(k1+3);
    end

    %% 求解
    T=A\b;
   
    %% 计算L2 误差
    % 计算 un 在求积节点处取值 
    C=zeros(1,1001);
    for k=1:1:N
        C=C+T(k)*(1-x).*x.^k;
    end
    E=(C-u).^2;

    % 利用复化Simpson公式计算L^2 误差
    err2(i)=sqrt(ComSimpson(0, 1, E));

    if(i==10)
       figure
       spy(A);
    end

    %% 计算条件数，cond(A, 2) 计算 A的从属 2 范数
    condB(i)=cond(A, 2 );
end

%% 画误差图
% loglog() 更适用于画指数函数来显示其变化曲线
figure
loglog(err1, '-.r', 'DisplayName', 'Basis I');
hold on
loglog(err2, '--k', 'DisplayName', 'Basis II');
legend('Show')
hold off


%% plot the condition number
figure
loglog(condA, '-.r' , 'DisplayName', 'Basis I');
hold on
loglog(condB, '--k', 'DisplayName', 'Basis II');
legend('Show')
hold off