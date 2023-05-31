clc
clear 
tic
%%创建PDE模型
model=createpde;
R1=[3,4,0,1,1,0,0,0,1,1]';
gn=[R1];
g=decsg(gn);
geometryFromEdges(model,g);

% 进行循环的次数
MK=10;

%%声明
err2=zeros(MK, 1);    errf=zeros(MK,1);
errH1=zeros(MK, 1);   errL2=zeros(MK,1);
CondN=zeros(MK, 1);
h0=[0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.015, 0.012, 0.01];

for ni=1:1:MK  

    %% 网格信息
    mesh=generateMesh(model,'GeometricOrder','linear','Hmax',h0(i));
    [p, e, t ]= meshToPet(mesh);

    % figure; 
    % pdemesh(p, e, t );

    Np=length(p);
    A=zeros(Np,Np);
    b=zeros(Np,1);
    u=zeros(Np,1);
    K=zeros(3,3);
    L=zeros(3);


    %% 给出变分问题 a(u, v)
    % 对于此问题 p(x)=1, q(x)=pi^2/4
    %p1=-2*pi^2; q1=1; q2=1;

    %% 定义标准区域上的数值求积公式，定义求积节点
    % 这里利用3点高斯公式
    xi=[2/3, 1/6, 1/6;
        1/6, 2/3, 1/6;
        1/6, 1/6, 2/3];
    w=[1/3, 1/3, 1/3];
    Nxi=length(xi);

    %% 组装刚度矩阵
    VecN=xi;
    VecNx=ones(Nxi,3);
    VecNy=ones(Nxi,3);

    for i=1:1:length(t)
        x1=p(1, t(1 ,i)); y1=p(2, t(1,i));
        x2=p(1, t(2 ,i)); y2=p(2, t(2,i));
        x3=p(1, t(3 ,i)); y3=p(2, t(3,i));
        S=det([ones(1,3);[x1, x2, x3]; [y1, y2, y3]])/2;
        DVecx=[y2-y3, y3-y1, y1-y2];
        DVecy=-[x3-x3, x3-x1; x1-x2];

        %求积节点
        tempx=xi*[x1 x2 x3]';
        tempy=xi*[y1 y2 y3]';

        %定义参数矩阵 p q1 q2
        VecP=2*ones(Nxi, 1)*pi^2;
        VecQ1=-ones(Nxi, 1);
        VecQ2=-ones(Nxi, 1);

        %定义右端函数 F
        VecF1=2*pi^2*tempx.*tempy;

        %形成刚度矩阵
        for j1=1:1:3
            for j2=1:1:3
                temp=VecP.*VecN(:, j1).*VecN(:, j2)...
                    +VecQ1.*VecNx(:, j1).*VecNx(:, j2)*DVecx(j1)*DVecx(j2)/(4*S^2)...
                    +VecQ2.*VecNy(:, j1).*VecNy(:, j2)*DVecy(j1)*DVecy(j2)/(4*S^2);
                K(j1,j2)=S*w*temp;

                A(t(j1,i),t(j2,i))=A(t(j1,i),t(j2,i))+K(j1,j2);
            end

            temp=VecF1.*VecN(:,j1);
            L(j1)=S*w*temp;

            b(t(j1,i))=b(t(j1,i))+L(j1);
        end
    end    

    %自然边界条件处理
    xii=GaussianQPoints(0, 1, 4);
    Nxii=[1-xii;xii];

    for i=1:1:length(e)
        x1=p(1, e(1, i));    y1=p(2, e(1, i));
        x2=p(1, e(2, i));    y2=p(2, e(2, i));
        hxy=sqrt((x1-x2)^2+(y1-y2)^2);          %边界处的 hxy=hx 或者 hy
        tempx=x1+(x2-x1)*xii;
        tempy=y1+(y2-y1)*xii;

        if(e(5,i)==2)
            VecF=-(tempy-pi*sin(pi*tempy));
            L=zeros(2,1);
            for j=1:1:2
                temp=VecF.*Nxii(j, :)*hxy;
                L(j)=GaussianQuadrature(0, 1, temp, 4);
                b(e(j,i))=b(e(j,i))+L(j);
            end
        end

        if(e(5,i)==3)
            VecF=-(tempx-pi*sin(pi*tempx));
            L=zeros(2,1);
            for j=1:1:2
                temp=VecF.*Nxii(j, :)*hxy;
                L(j)=GaussianQuadrature(0, 1, temp, 4);
                b(e(j,i))=b(e(j,i))+L(j);
            end
        end      
    end
    
    % 处理本质边界条件
    for i=1:1:length(e)

        if(e(5,i)==1)
            for j=1:1:Np
                A(e(1,i),j)=0;
            end
            A(e(1,i),e(1,i))=1;
        end

        if(e(5,i)==4)
            for j=1:1:Np
                A(e(1,i),j)=0;
            end
            A(e(1,i),e(1,i))=1;
        end

    end

    A=sparse(A);

    %% 求解系数矩阵
    u=A\b;
    CondN(ni)=condest(A);

    %% 计算 L2 和 H1 误差
    % 利用复化Guass公式分区间段进行数值积分
    % 每个区间段内 u=u_{i1} N0(\xi)+ u_{ic} Nc(\xi)  + u_{i2} N1(\xi)
    % 每个区间段内 u'=u_{i1} N0'(\xi)+ u_{ic} Nc'(\xi)  + u_{i2} N1'(\xi)
    for i=1:1:length(t)
        
        x1=p(1, t(1 ,i)); y1=p(2, t(1,i));
        x2=p(1, t(2 ,i)); y2=p(2, t(2,i));
        x3=p(1, t(3 ,i)); y3=p(2, t(3,i));
        S=det([ones(1,3);[x1, x2, x3]; [y1, y2, y3]])/2;
        DVecx=[y2-y3, y3-y1, y1-y2];
        DVecy=-[x3-x3, x3-x1; x1-x2];

        %求积节点
        tempx=xi*[x1 x2 x3]';
        tempy=xi*[y1 y2 y3]';
   

        % u 和 真解的值
        tempu=[u(t(1 ,i)) u(t(2 ,i)) u(t(3 ,i))]*xi;
        tempv=tempx.*tempy+sin(pi*tempx).*sin(pi*tempy);

        % u 和 真解的导数值
        tempdux=[u(t(1 ,i)) u(t(2 ,i)) u(t(3 ,i))]*[DVecy' DVecy' DVecy']/(2*S);
        tempduy=[u(t(1 ,i)) u(t(2 ,i)) u(t(3 ,i))]*[DVecx' DVecx' DVecx']/(2*S);
        tempdvx=tempy+pi*cos(pi*tempx).*sin(pi*tempy);
        tempdvy=tempx+pi*sin(pi*tempx).*cos(pi*tempy);

        % 计算 L2 误差
        temp=(tempu-tempv).^2*S;
        err2(ni)=err2(ni)+GaussianQuadrature(0, 1, temp, 4);

        % 计算 H1 半范数
        temp=((tempdux-tempdvx).^2 + (tempduy-tempdvy).^2)*S;
        errf(ni)=errf(ni)+GaussianQuadrature(0, 1, temp, 4);


    end

    % 计算 L2 和 H1 范数
    errH1(ni)=sqrt(errf(ni)+err2(ni));
    errL2(ni)=sqrt(err2(ni));
end

%% 画数值解
v=zeros(Np)
X=p(1,:);
Y=p(2,:);
v=X.*Y+sin(pi*X).*sin(pi*Y);

figure
plot3(X,Y,u);
hold on;
plot3(X,Y,v);
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

