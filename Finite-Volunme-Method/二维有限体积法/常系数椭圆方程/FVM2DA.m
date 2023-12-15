%% 二维矩形剖分上的差分法：逐点形成（有限体积法）

%% 测试问题 u*=xy+sin(pi x) sin(pi y)
% \Delta u+2pi^2 u=2pi^2 xy in [0, 1]*[0, 1]
% u(0, y)=u(x, 0)=0
% \partial_x u=y-\pi sin \pi y \in {x=1}*[0, 1]
% \partial_y u=x-\pi sin \pi x \in [0, 1]*{y=1}

clear

%% 区域 \Omega=[a1, a2]*[a3, a4]
a1=0; a2=1; a3=0; a4=1;

NK=15;
errC=zeros(NK, 1); errL=zeros(NK, 1); H0err=zeros(NK, 1); meshsize=zeros(NK, 1);

for ni=1:1:NK
    N=5*ni; M=10*ni;
    % 声明总刚度矩阵
    A=zeros((N+1)*(M+1), (N+1)*(M+1));
    b=zeros((N+1)*(M+1), 1);
    
    % 依旧考虑均匀网格
    hx=(a2-a1)/M;  hy=(a4-a3)/N;
    
    % 界点 Gamma1 \cup Gamma4
    Gamma1=[1:1:M+1];
    Gamma4=1+[0:1:N]*(M+1);
    
    % 界点 Gamma2 \cup Gamma3
    Gamma2=(M+1)*[1:1:N+1];
    Gamma3=N*(M+1)+[1:1:M+1];    
    
    %% Step 3: 形成系数矩阵
    % 内点矩阵
    for i=2:1:M
        for j=2:1:N
            k=(j-1)*(M+1)+i;
            tempx=a1+(i-1)*hx;
            tempy=a3+(j-1)*hy;
    
            A(k, k)=-2/hy^2-2/hx^2+2*pi^2;
            A(k, k-1)=1/hx^2;
            A(k, k+1)=1/hx^2;
            A(k, k+M+1)=1/hy^2;
            A(k, k-M-1)=1/hy^2;
            b(k)=2*pi^2*tempx*tempy;
        end
    end
    
    % 界点的处理:    Gamma 2
    for i=2:1:length(Gamma2)-1
        tempx=a2;
        tempy=(i-1)*hy+a3;
        tempy1=tempy+0.5*hy;
        tempy2=tempy-0.5*hy;
    
        A(Gamma2(i), Gamma2(i))=-1/hy^2-1/hx^2+2*pi^2/2;
        A(Gamma2(i), Gamma2(i)-1)=1/hx^2;
        A(Gamma2(i), Gamma2(i)+M+1)=1/(2*hy^2);
        A(Gamma2(i), Gamma2(i)-M-1)=1/(2*hy^2);
        b(Gamma2(i))=0.5*2*pi^2*tempx*tempy-1/hx*(tempy-pi*sin(pi*tempy));
    end
    
    % 界点的处理:    Gamma 3
    for i=2:1:length(Gamma3)-1
        tempx=(i-1)*hx+a1;
        tempy=a4;
        tempx1=tempx+0.5*hx;
        tempx2=tempx-0.5*hx;
    
        A(Gamma3(i), Gamma3(i))=-1/hy^2-1/hx^2+2*pi^2/2;
        A(Gamma3(i), Gamma3(i)-1)=1/(2*hx^2);
        A(Gamma3(i), Gamma3(i)+1)=1/(2*hx^2);
        A(Gamma3(i), Gamma3(i)-M-1)=1/hy^2;
        b(Gamma3(i))=0.5*2*pi^2*tempx*tempy-1/hy*(tempx-pi*sin(pi*tempx));
    end
    
    % 角点 Gamma2 \cap Gamma3
    tempx=a2;   tempy=a4;
    tempx2=a2-0.5*hx;
    tempy2=a4-0.5*hy;
    A((N+1)*(M+1), (N+1)*(M+1))=-1/(2*hx^2)-1/(2*hy^2)+2*pi^2/4;
    A((N+1)*(M+1), (N+1)*(M+1)-1)=1/(2*hx^2);
    A((N+1)*(M+1), (N+1)*(M+1)-M-1)=1/(2*hy^2);
    b((N+1)*(M+1))=0.25*2*pi^2*tempx*tempy-0.5/hx*(a4-pi*sin(pi*a4))-0.5/hy*(a2-pi*sin(pi*a2));
    
    % 界点处理：Gamma1
    A(Gamma1, :)=0;
    for i=1:1:length(Gamma1)
        A(Gamma1(i), Gamma1(i))=1;
        b(Gamma1(i))=0;
    end
    
    % 界点处理：Gamma4
    A(Gamma4, :)=0;
    for i=1:1:length(Gamma4)
        A(Gamma4(i), Gamma4(i))=1;
        b(Gamma4(i))=0;
    end
    
    %% Step 4: 系数矩阵求解
    u=A\b;
    
    %% Step 5: 误差计算
    utrue=zeros((N+1)*(M+1), 1);
    for i=1:1:M+1
      tempx=a1+(i-1)*hx;
      for j=1:1:N+1
          k=(j-1)*(M+1)+i;
          tempy=a3+(j-1)*hy;
          utrue(k)=tempx*tempy+sin(pi*tempx)*sin(pi*tempy);
      end
    end
    
    % C 范数
    errC(ni)=max(abs(utrue-u));
    
    %  L2 范数
    errL(ni)=(utrue-u)'*(utrue-u)*hx*hy;
    errL(ni)=sqrt(errL(ni));

    % 网格最大边长
    meshsize(ni)=sqrt(hx^2+hy^2);
end

%% 作图
loglog(meshsize, errC, '-.', 'DisplayName', 'C-norm error');
hold on
loglog(meshsize, meshsize.^2, '--k', 'DisplayName', 'O(h^2)')
legend('Show')
hold off

figure
loglog(meshsize, errL, '-.', 'DisplayName', '0-norm error');
hold on
loglog(meshsize, meshsize.^2, '--k', 'DisplayName', 'O(h^2)')
legend('Show')
hold off


% 数值阶计算
alpha1=zeros(NK-1, 1);      % C 范数
alpha2=zeros(NK-1, 1);      % L 范数
for i=1:1:NK-1
    alpha1(i)=(log(errC(i+1))-log(errC(i)))/(log(meshsize(i+1))-log(meshsize(i)));
    alpha2(i)=(log(errL(i+1))-log(errL(i)))/(log(meshsize(i+1))-log(meshsize(i)));
end

