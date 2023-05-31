clc
clear 
tic

% 进行循环的次数
MK=1;

errC=zeros(MK, 1);    err0=zeros(MK,1);     CondN=zeros(MK, 1);
f  =@(x,y)(2*pi^2*x*y);
q  =@(x,y)(2*pi^2);
dux=@(y)(y-pi*sin(pi*y));
duy=@(x)(x-pi*sin(pi*x));

for ni=1:1:MK
    n=5*ni;
    a=0; b=1; 
    y=zeros((n+1)^2, 1);
    u=zeros((n+1)^2, 1);
    p=zeros(n+1,1); 
    I=zeros(n, 3);
    ix=ones(5*(n+1)^2,1);%稀疏矩阵行指标
    jy=ones(5*(n+1)^2,1);%        列指标
    vz=zeros(5*(n+1)^2,1);%       函数值
    vz(1)=1;
    
    %% 网格信息
    % 这里考虑均匀剖分 
    h=(b-a)/n;
    p(:, 1)=a+(0:1:n)*h;
    I(:, 1)=(1:1:n); I(:, 2)=(2:1:n+1); I(:, 3)=(3:1:n+2);


    %% 计算系数矩阵
    for i=1:1:n-1
        for j=1:1:n-1
            x1=p(I(i,1)); x2=p(I(i,2)); x3=p(I(i,3));
            y1=p(I(j,1)); y2=p(I(j,2)); y3=p(I(j,3));
            hx=x2-x1;     hy=y2-y1;

            k =(I(i,2)-1)*(n+1)+I(j,2);
            m1=(I(i,3)-1)*(n+1)+I(j,2);
            m2=(I(i,1)-1)*(n+1)+I(j,2);
            m3=(I(i,2)-1)*(n+1)+I(j,3);
            m4=(I(i,2)-1)*(n+1)+I(j,1);           

            i1=5*(i*(n+1)+j);
            ix(i1+(1:1:5))=k;
            jy(i1+1)=m1; jy(i1+2)=m2; jy(i1+3)=k;  jy(i1+4)=m3;  jy(i1+5)=m4;
            vz(i1+1)=+1/hx^2;
            vz(i1+2)=+1/hx^2;
            vz(i1+3)=-2/hx^2-2/hy^2+q(x2,y2);
            vz(i1+4)=+1/hy^2;
            vz(i1+5)=+1/hy^2;

            y(k)=y(k)+ f(x2,y2);
            
        end
    end
       
    % 处理边值条件
    i=n;
    for j=1:1:n-1
        x1=p(I(i,1)); x2=p(I(i,2)); 
        y1=p(I(j,1)); y2=p(I(j,2)); y3=p(I(j,3));
        hx=x2-x1;     hy=y2-y1;

        k =(I(i,2)-1)*(n+1)+I(j,2);
        m1=(I(i,3)-1)*(n+1)+I(j,2);
        m2=(I(i,1)-1)*(n+1)+I(j,2);
        m3=(I(i,2)-1)*(n+1)+I(j,3);
        m4=(I(i,2)-1)*(n+1)+I(j,1);
        
        i1=5*(i*(n+1)+j);
        ix(i1+(1:1:5))=k;
        jy(i1+1)=m1; jy(i1+2)=m2; jy(i1+3)=k;  jy(i1+4)=m3;  jy(i1+5)=m4;
  
        vz(i1+2)= +1/hx^2;
        vz(i1+3)= -1/hx^2-1/hy^2+q(x2,y2)/2;
        vz(i1+4)= +1/hy^2/2;
        vz(i1+5)= +1/hy^2/2;

        y(k) = y(k)+ f(x2,y2)/2- dux(y2)/hx;
    end 

    j=n;
    for i=1:1:n-1
        x1=p(I(i,1)); x2=p(I(i,2)); x3=p(I(i,3));
        y1=p(I(j,1)); y2=p(I(j,2)); 
        hx=x2-x1;     hy=y2-y1;

        k =(I(i,2)-1)*(n+1)+I(j,2);
        m1=(I(i,3)-1)*(n+1)+I(j,2);
        m2=(I(i,1)-1)*(n+1)+I(j,2);
        m3=(I(i,2)-1)*(n+1)+I(j,3);
        m4=(I(i,2)-1)*(n+1)+I(j,1);

        i1=5*(i*(n+1)+j);
        ix(i1+(1:1:5))=k;
        jy(i1+1)=m1; jy(i1+2)=m2; jy(i1+3)=k;  jy(i1+4)=m3;  jy(i1+5)=m4;
        
        vz(i1+1)= +1/hx^2/2;
        vz(i1+2)= +1/hx^2/2;
        vz(i1+3)= -1/hy^2-1/hx^2+q(x2,y2)/2;
        vz(i1+5)= +1/hy^2;

        y(k) = f(x2,y2)/2- duy(x2)/hy;

    end

    i=n; j=n;
    x1=p(I(i,1)); x2=p(I(i,2)); 
    y1=p(I(j,1)); y2=p(I(j,2));
    hx=x2-x1;     hy=y2-y1;

    k =(I(i,2)-1)*(n+1)+I(j,2);
    m1=(I(i,3)-1)*(n+1)+I(j,2);
    m2=(I(i,1)-1)*(n+1)+I(j,2);
    m3=(I(i,2)-1)*(n+1)+I(j,3);
    m4=(I(i,2)-1)*(n+1)+I(j,1);
    i1=5*(i*(n+1)+j);
    ix(i1+(1:1:5))=k;
    jy(i1+1)=m1; jy(i1+2)=m2; jy(i1+3)=k;  jy(i1+4)=m3;  jy(i1+5)=m4;

    y(k) = f(x2,y2)/4- duy(x2)/hy/2- dux(y2)/hx/2;
    vz(i1+2)= +1/hx^2/2;
    vz(i1+3)= -1/hy^2/2-1/hx^2/2+q(x2,y2)/4;
    vz(i1+5)= +1/hy^2/2;

    y(1:1:n+1)=0;
    y((1:1:n)*(n+1)+1)=0;
    
    
    j=0;
    for i=0:1:n
        % A(i,:)=0  A(i,i)=1;
        i1=5*(i*(n+1)+j);
        ix(i1+3)=i+1;
        jy(i1+3)=i+1;
        vz(i1+3)=1;
           
    end

    i=0;
    for j=0:1:n
        %A((i-1)*(n+1)+1,(i-1)*(n+1)+1)=1;
        i1=5*(i*(n+1)+j);
        ix(i1+3)=j*(n+1)+1;
        jy(i1+3)=j*(n+1)+1;
        vz(i1+3)=1;
    end

    A=sparse(ix,jy,vz,(n+1)^2,(n+1)^2);
    A=A(:,(3:(n+1)^2+2));
    %% 求解系数矩阵
    
    u=A\y;

    u(1:1:n+1)=0;
    u((0:1:n)*(n+1)+1)=0;

    CondN(ni)=condest(A);

    %% 计算误差
    v=zeros((n+1)^2,1);
    for i=1:1:n+1
        for j=1:1:n+1
            v((i-1)*(n+1)+j)=(i-1)*(j-1)/n/n+sin(pi*(i-1)/n)*sin(pi*(j-1)/n);
        end
    end

    temp=u-v;
    errC(ni)=max(abs(temp));
    err0(ni)=sqrt(sum(temp.^2))*h;

end

toc
%% 画数值解

U=zeros(n+1,n+1);
V=zeros(n+1,n+1);
[X,Y]=meshgrid(0:1/n:1);
V=X.*Y+sin(pi*X).*sin(pi*Y);
for i=1:1:n+1
    for j=1:1:n+1
        U(i,j)=u((i-1)*(n+1)+j);
    end
end
V=X.*Y+sin(pi*X).*sin(pi*Y);

figure
mesh(X,Y,U);
hold on;
mesh(X,Y,V);
hold off


%% 画误差图
% C 误差
N=5*(1:1:MK);
h=1./N;
figure
loglog(h, errC, '--r', 'DisplayName', 'errC');
hold on
loglog(h,h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

% H0 误差
figure
loglog(h, err0, '--r', 'DisplayName', 'err0');
hold on
loglog(h, h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

%% Condition number
figure
loglog(h, CondN, '--r', 'DisplayName', 'Condition Number');
hold on
loglog(h, 1/4*h.^-3, '-.k', 'DisplayName', 'y=1/h^{-3}');
legend('Show');
hold off

% 数值阶计算
alpha1=zeros(MK-1, 1);
alpha2=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errC(i+1))-log(errC(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(err0(i+1))-log(err0(i)))/(log(h(i+1))-log(h(i)));    
end

c=ones(1e4,1);
a=eye(1e5);
tic
b=a\c;
toc
tic
a=sparse(a);
b=a\c;
toc

