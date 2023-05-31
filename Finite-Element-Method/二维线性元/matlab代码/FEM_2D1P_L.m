clc
clear 
tic

% 进行循环的次数
MK=10;

err2=zeros(MK, 1);    errf=zeros(MK,1);
errH1=zeros(MK, 1);   errL2=zeros(MK,1);
CondN=zeros(MK, 1);

for ni=1:1:MK
    n=5*ni;
    a=0; b=1;
    A=zeros((n+1)^2, (n+1)^2); 
    y=zeros((n+1)^2, 1);
    u=zeros((n+1)^2, 1);
    p=zeros(1,n+1); 
    I=zeros(2, n);
    
    
    
    %% 网格信息
    % 这里考虑均匀剖分 
    h=(b-a)/n;
    p(1, : )=a+(0:1:n)*h;
    I(1, : )=(1:1:n); I(2, : )=(2:1:n+1);


    %% 给出变分问题 a(u, v)
    % 对于此问题 p(x)=1, q(x)=pi^2/4
    p1=-2*pi^2; q1=1; q2=1;

    %% 定义标准区域上的数值求积公式，定义求积节点
    % 这里利用复化ComSimson公式

    M=8; hm=1/(2*M); 
    xi=zeros(2*M+1,1);  xi=(0:1:2*M)'*hm;
    yi=zeros(2*M+1,1);  yi=(0:1:2*M)'*hm;
%     M=1;
%     xi=zeros(3,1); 
%     xi(1)=(b+a-sqrt(3/5)*(b-a))/2;
%     xi(2)=(b+a)/2;
%     xi(3)=(b+a+sqrt(3/5)*(b-a))/2;
% 
%     yi=zeros(3,1); 
%     yi(1)=(b+a-sqrt(3/5)*(b-a))/2;
%     yi(2)=(b+a)/2;
%     yi(3)=(b+a+sqrt(3/5)*(b-a))/2;

    II=ones(2*M+1,1);
    temp=zeros(2*M+1,2*M+1);

    %% 计算I_{i，j)单元刚度矩阵
    for i=1:1:n
        for j=1:1:n
            i1=I(1, i); i2=I(2, i); x1=p(i1);  x2=p(i2); hi=x2-x1;
            j1=I(1, j); j2=I(2, j); y1=p(j1);  y2=p(j2); hj=y2-y1;
            m1=(i1-1)*(n+1)+j1;
            m2=(i1-1)*(n+1)+j2;
            m3=(i2-1)*(n+1)+j1;
            m4=(i2-1)*(n+1)+j2;
            
            
            % 计算
            for k=1:1:2*M+1
                temp(k, : )=hi*hj*p1*(1-xi(k))^2*(1-yi).^2+q1*(1-yi).^2+q2*(1-xi(k))^2;                       
            end
            a11=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=hi*hj*p1*(1-xi(k))^2*(1-yi).*yi+q1*yi.*(1-yi)-q2*(1-xi(k))^2;                       
            end     
            a12=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=hi*hj*p1*(1-xi(k))*xi(k)*(1-yi).^2-q1*(1-yi).^2+q2*(1-xi(k))*xi(k);                       
            end       
            a13=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=hi*hj*p1*(1-xi(k))*xi(k)*(1-yi).*yi-q1*yi.*(1-yi)-q2*(1-xi(k))*xi(k);                       
            end     
            a14=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=hi*hj*p1*(1-xi(k))^2*yi.^2+q1*yi.^2+q2*(1-xi(k))^2;                       
            end                        
            a22=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=hi*hj*p1*(1-xi(k))*xi(k)*yi.*(1-yi)-q1*yi.*(1-yi)-q2*(1-xi(k))*xi(k);                       
            end          
            a23=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=hi*hj*p1*(1-xi(k))*xi(k)*yi.^2-q1*yi.^2+q2*(1-xi(k))*xi(k);                       
            end                  
            a24=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=hi*hj*p1*xi(k)^2*(1-yi).^2 + q1*(1-yi).^2 + q2*xi(k)^2;                       
            end                          
            a33=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=hi*hj*p1*xi(k)^2*(1-yi).*yi+q1*yi.*(1-yi)-q2*xi(k)^2;                       
            end  
            a34=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=hi*hj*p1*xi(k)^2*yi.^2+q1*yi.^2+q2*xi(k)^2;                       
            end                              
            a44=ComSimpson(temp);

            a21=a12;
            a31=a13;    a32=a23;
            a41=a14;    a42=a24;    a43=a34;
            
            
            % 计算 b1, b2, b3, b4
            for k=1:1:2*M+1
                temp(k, : )=-hi*hj*2*pi^2*(x1+hi*xi(k))*(y1+hj*yi)*(1-xi(k)).*(1-yi);                       
            end
            b1=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=-hi*hj*2*pi^2*(x1+hi*xi(k))*(y1+hj*yi)*(1-xi(k)).*yi;                       
            end
            b2=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=-hi*hj*2*pi^2*(x1+hi*xi(k))*(y1+hj*yi)*xi(k).*(1-yi);                       
            end
            b3=ComSimpson(temp);


            for k=1:1:2*M+1
                temp(k, : )=-hi*hj*2*pi^2*(x1+hi*xi(k))*(y1+hj*yi)*xi(k).*yi;                       
            end
            b4=ComSimpson(temp);
            
            

            % 组装刚度矩阵
            A(m1, m1)=A(m1, m1)+a11;  A(m1, m2)=A(m1, m2)+a12;  A(m1, m3)=A(m1, m3)+a13;  A(m1, m4)=A(m1, m4)+a14;
            A(m2, m1)=A(m2, m1)+a21;  A(m2, m2)=A(m2, m2)+a22;  A(m2, m3)=A(m2, m3)+a23;  A(m2, m4)=A(m2, m4)+a24;
            A(m3, m1)=A(m3, m1)+a31;  A(m3, m2)=A(m3, m2)+a32;  A(m3, m3)=A(m3, m3)+a33;  A(m3, m4)=A(m3, m4)+a34;
            A(m4, m1)=A(m4, m1)+a41;  A(m4, m2)=A(m4, m2)+a42;  A(m4, m3)=A(m4, m3)+a43;  A(m4, m4)=A(m4, m4)+a44;

            y(m1, 1)=y(m1, 1)+b1;     y(m2, 1)=y(m2, 1)+b2;     y(m3, 1)=y(m3, 1)+b3;     y(m4, 1)=y(m4, 1)+b4;
        end
    end
    
    
    % 处理本质边界条件
    i=n;
    for j=1:1:n
        i1=I(1, i); i2=I(2, i); x1=p(i1);  x2=p(i2); hi=x2-x1;
        j1=I(1, j); j2=I(2, j); y1=p(j1);  y2=p(j2); hj=y2-y1;
        m1=(i1-1)*(n+1)+j1;
        m2=(i1-1)*(n+1)+j2;
        m3=(i2-1)*(n+1)+j1;
        m4=(i2-1)*(n+1)+j2;

        b5=hj*(y1+hj*yi-pi*sin(pi*(y1+hj*yi))).*(1-yi);
        b6=hj*(y1+hj*yi-pi*sin(pi*(y1+hj*yi))).*yi;

        b=ComSimpsonX(b5);
        y(m3, 1)=y(m3, 1)+b;

        b=ComSimpsonX(b6);
        y(m4, 1)=y(m4, 1)+b;


    end

    j=n;
    for i=1:1:n
        i1=I(1, i); i2=I(2, i); x1=p(i1);  x2=p(i2); hi=x2-x1;
        j1=I(1, j); j2=I(2, j); y1=p(j1);  y2=p(j2); hj=y2-y1;
        m1=(i1-1)*(n+1)+j1;
        m2=(i1-1)*(n+1)+j2;
        m3=(i2-1)*(n+1)+j1;
        m4=(i2-1)*(n+1)+j2;

        b5=hi*(x1+hi*xi-pi*sin(pi*(x1+hi*xi))).*(1-xi);
        b6=hi*(x1+hi*xi-pi*sin(pi*(x1+hi*xi))).*xi;

        b=ComSimpsonX(b5);
        y(m2, 1)=y(m2, 1)+b;

        b=ComSimpsonX(b6);
        y(m4, 1)=y(m4, 1)+b;

    end


    y(1:1:n+1)=0;
    y((0:1:n)*(n+1)+1)=0;
    for i=1:1:n+1
        A(i,:)=0;
        A(i,i)=1;
        A((i-1)*(n+1)+1,:)=0;
        A((i-1)*(n+1)+1,(i-1)*(n+1)+1)=1;
    end

    A=sparse(A);
    %% 求解系数矩阵
    
    u=A\y;

    u(1:1:n+1)=0;
    u((0:1:n)*(n+1)+1)=0;

    CondN(ni)=condest(A);

    %% 计算 L2 和 H1 误差
    % 利用复化Guass公式分区间段进行数值积分
    % 每个区间段内 u=u_{i1} N0(\xi)+ u_{ic} Nc(\xi)  + u_{i2} N1(\xi)
    % 每个区间段内 u'=u_{i1} N0'(\xi)+ u_{ic} Nc'(\xi)  + u_{i2} N1'(\xi)
    for i=1:1:n
        for j=1:1:n
            i1=I(1, i); i2=I(2, i); x1=p(i1);  x2=p(i2); hi=x2-x1;
            j1=I(1, j); j2=I(2, j); y1=p(j1);  y2=p(j2); hj=y2-y1;
            m1=(i1-1)*(n+1)+j1;
            m2=(i1-1)*(n+1)+j2;
            m3=(i2-1)*(n+1)+j1;
            m4=(i2-1)*(n+1)+j2;

            % u 和 真解的值
            tempu=u(m1)*(1-xi)*(1-yi)'+ u(m2)*(1-xi)*yi' + u(m3)*xi*(1-yi)' + u(m4)*xi*yi';
            tempv=(x1+hi*xi)*(y1+hj*yi)'+sin(pi*(x1+hi*xi))*sin(pi*(y1+hj*yi))';

            % u 和 真解的导数值
            tempdux=-u(m1)*II*(1-yi)'/hi -u(m2)*II*yi'/hi +u(m3)*II*(1-yi)'/hi +u(m4)*II*yi'/hi;
            tempduy=-u(m1)*(1-xi)*II'/hj +u(m2)*(1-xi)*II'/hj -u(m3)*xi*II'/hj +u(m4)*xi*II'/hj;
            tempdvx=II*(y1+hj*yi)'+pi*cos(pi*(x1+hi*xi))*sin(pi*(y1+hj*yi))';
            tempdvy=(x1+hi*xi)*II'+pi*sin(pi*(x1+hi*xi))*cos(pi*(y1+hj*yi))';

            % 计算 L2 误差
            temp=(tempu-tempv).^2*hi*hj;
            err2(ni)=err2(ni)+ComSimpson(temp);

            % 计算 H1 半范数
            temp=((tempdux-tempdvx).^2 + (tempduy-tempdvy).^2)*hi*hj;
            errf(ni)=errf(ni)+ComSimpson(temp);

        end
    end

    % 计算 L2 和 H1 范数
    errH1(ni)=sqrt(errf(ni)+err2(ni));
    errL2(ni)=sqrt(err2(ni));
end

%% 画数值解

U=zeros(n+1,n+1);
V=zeros(n+1,n+1);
[X,Y]=meshgrid([0:1/n:1]);
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
% L^2 误差
N=5*[1:1:MK];
h=1./N;
figure
loglog(h, errL2, '--r', 'DisplayName', 'L^2 error');
hold on
loglog(h,h.^2/4, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

% H1 误差
figure
loglog(h, errH1, '--r', 'DisplayName', 'H1 error');
hold on
loglog(h, h*10, '-.k', 'DisplayName', 'O(h)')
legend('Show');
hold off

%% Condition number
figure
loglog(h, CondN, '--r', 'DisplayName', 'Condition Number');
hold on
loglog(h, 1/4*h.^-2, '-.k', 'DisplayName', 'y=1/h^{-2}');
legend('Show');
hold off

% 数值阶计算
alpha1=zeros(MK-1, 1);
alpha2=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errL2(i+1))-log(errL2(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(errH1(i+1))-log(errH1(i)))/(log(h(i+1))-log(h(i)));    
end
toc

