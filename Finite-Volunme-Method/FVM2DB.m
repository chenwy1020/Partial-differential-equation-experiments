%% 二维矩形剖分上的差分法：单元刚度矩阵方式

%% 测试问题 u*=xy+sin(pi x) sin(pi y)
% \Delta u+2pi^2 u=2pi^2 xy in [0, 1]*[0, 1]
% u(0, y)=u(x, 0)=0
% \partial_x u=y-\pi sin \pi y \in {x=1}*[0, 1]
% \partial_y u=x-\pi sin \pi x \in [0, 1]*{y=1}

clear

%% 区域 \Omega=[a1, a2]*[a3, a4]
a1=0; a2=1; a3=0; a4=1;

NK=15;
errC=zeros(NK, 1); errL=zeros(NK, 1); meshsize=zeros(NK, 1);

for ni=1:1:NK
    N=5*ni; M=10*ni;
    % 声明总刚度矩阵
    A=zeros((N+1)*(M+1), (N+1)*(M+1));
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

    %% Step 3: 生成系数矩阵
    for i=1:1:N*M
        i1=I(1, i); i2=I(2, i); i3=I(3, i); i4=I(4, i);

        A(i1, i1)=A(i1, i1)-1/(2*hy^2)-1/(2*hx^2)+0.25*2*pi^2;
        A(i1, i2)=A(i1, i2)+1/(2*hx^2);
        A(i1, i3)=A(i1, i3)+1/(2*hy^2);
        b(i1)=b(i1)+0.25*2*pi^2*p(1, i1)*p(2, i1);

        A(i2, i2)=A(i2, i2)-1/(2*hy^2)-1/(2*hx^2)+0.25*2*pi^2;
        A(i2, i1)=A(i2, i1)+1/(2*hx^2);
        A(i2, i4)=A(i2, i4)+1/(2*hy^2);
        b(i2)=b(i2)+0.25*2*pi^2*p(1, i2)*p(2, i2);

        A(i3, i3)=A(i3, i3)-1/(2*hy^2)-1/(2*hx^2)+0.25*2*pi^2;
        A(i3, i1)=A(i3, i1)+1/(2*hy^2);
        A(i3, i4)=A(i3, i4)+1/(2*hx^2);
        b(i3)=b(i3)+0.25*2*pi^2*p(1, i3)*p(2, i3);

        A(i4, i4)=A(i4, i4)-1/(2*hy^2)-1/(2*hx^2)+0.25*2*pi^2;
        A(i4, i2)=A(i4, i2)+1/(2*hy^2);
        A(i4, i3)=A(i4, i3)+1/(2*hx^2);
        b(i4)=b(i4)+0.25*2*pi^2*p(1, i4)*p(2, i4);
    end

    % 自然边界条件的处理
    for i=1:1:2*(M+N)
        i1=E(1, i); i2=E(2, i);
        if(E(3, i)==2)
            tempy1=p(2, i1);
            tempy2=p(2, i2);
            tempy=0.5*(p(2, i1)+p(2, i2));
            b(i1)=b(i1)-(tempy1-pi*sin(pi*tempy1))/hx/2;
            b(i2)=b(i2)-(tempy2-pi*sin(pi*tempy2))/hx/2;           
        elseif (E(3, i)==3)
            tempx1=p(1, i1);
            tempx2=p(1, i2);
            tempx=0.5*(p(1, i1)+p(1, i2));
            b(i1)=b(i1)-(tempx1-pi*sin(pi*tempx1))/hy/2;
            b(i2)=b(i2)-(tempx2-pi*sin(pi*tempx2))/hy/2;       
        end
    end

    % 本质边界条件的处理
    for i=1:1:2*(M+N)
        i1=E(1, i); i2=E(2, i);
        if(E(3, i)==1 || E(3, i)==4)
            A(i1, :)=0; A(i1, i1)=1;    b(i1)=0;
            A(i2, :)=0; A(i2, i2)=1;    b(i2)=0;
        end
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


