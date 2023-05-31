%% 二维三角网格剖分Lagrange 1次元

%% 测试问题 u*=xy+sin(pi x) sin(pi y)
% \Delta u+2pi^2 u=xy in [0, 1]*[0, 1]
% u(0, y)=u(x, 0)=0
% \partial_x u=y-\pi sin \pi y \in {x=1}*[0, 1]
% \partial_y u=x-\pi sin \pi x \in [0, 1]*{y=1}

clear

%% 网格生成: 构造 geometry description matrix
%   具体详见例：GenerateMeshA.m
model=createpde;
R1=[3, 4, 0, 1, 1, 0, 0 , 0 , 1, 1]';
gm=[R1];
g=decsg(gm);
geometryFromEdges(model,g); 

%% 声明变量
NK=1;
L2err=zeros(NK, 1); H1err=zeros(NK, 1); H0err=zeros(NK, 1);
AveH=zeros(NK, 1);
Aveh0=[0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.015, 0.012, 0.01]; 

for ni=1:1:NK

    %% Step2: 生成网格
    %   具体详见例：GenerateMeshA.m    
    mesh=generateMesh(model, 'GeometricOrder', 'linear', 'Hmax', Aveh0(ni));
    [p, e, t]=meshToPet(mesh);
    figure; pdemesh(p, e, t);

    G=t;
    
    % 声明总刚度矩阵
    Np=length(p);
    A=zeros(Np,Np );
    b=zeros(Np, 1);
    
    %% Step 3:  Lagrange 型元
    % 在每个三角单元内，节点按照逆时针排列
    
    % 高斯型求积公式
    [ w, xi ] = GaussianQuadraturePW(6);
    Nxi=length(w);
    
    VecN=xi;
    VecNx=ones(Nxi, 3);
    VecNy=ones(Nxi, 3);
    
    %% Step 4: 形成刚度矩阵
    % 此处采用不区分本质边界条件的方式
    % a(u, v)=\int_{\Omega} -q1ux vx-q2uy vy+puv...
    %   =int_{Omega} f1v-\int_{Gamma2}f2 vdy-\int_{Gamma3} f3 vdx
    for i=1:1:length(t)
        % Step 4.2 计算单元刚度矩阵
        K=zeros(3,3);   L=zeros(3, 1);
        x1=p(1, G(1, i));   y1=p(2, G(1, i));
        x2=p(1, G(2, i));   y2=p(2, G(2, i));
        x3=p(1, G(3, i));   y3=p(2, G(3, i));
        S=det([ones(1, 3); [x1, x2, x3]; [y1, y2, y3]])/2;
        DVecx=[y2-y3, y3-y1, y1-y2];
        DVecy=-[x2-x3, x3-x1, x1-x2];
    
        % 求积节点
        tempx=xi*[x1;x2;x3];    tempy=xi*[y1;y2;y3];
    
        % 定义参数矩阵 p, q1, q2
        VecP=2*pi^2*ones(Nxi,1);
        VecQ1=-ones(Nxi, 1);
        VecQ2=-ones(Nxi, 1);
    
        % 定义右端函数 F
        VecF1=2*pi^2*tempx.*tempy;
    
        % 形成单元刚度矩阵
        for j1=1:1:3
            for j2=1:1:3
                temp=VecP.*VecN(:, j1).*VecN(:, j2)...
                        +VecQ1.*VecNx(:, j1).*VecNx(:, j2)*DVecx(j1)*DVecx(j2)/(4*S^2)...
                        +VecQ2.*VecNy(:, j1).*VecNy(:, j2)*DVecy(j1)*DVecy(j2)/(4*S^2);
                K(j1, j2)=S*w*temp;
    
            end
            temp=VecF1.*VecN(:, j1);
            L(j1)=S*w*temp;
        end
    
       %  Step 4.3 组装刚度矩阵
       for j1=1:1:3
           for j2=1:1:3
               A(G(j1, i), G(j2, i))=A(G(j1, i), G(j2, i))+K(j1, j2);
           end
           b(G(j1, i), 1)=b(G(j1, i), 1)+L(j1);
       end
    end
    
    % Step 4.3 自然边界条件的处理
    % 确定一维的求积节点
    eta = GaussianQPoints(0, 1, 4);
    NE=[1-eta; eta]; 
    
    for i=1:1:length(e)
        x1=p(1, e(1, i));    y1=p(2, e(1, i));
        x2=p(1, e(2, i));    y2=p(2, e(2, i));
        hxy=sqrt((x1-x2)^2+(y1-y2)^2);
        tempx=x1+(x2-x1)*eta;
        tempy=y1+(y2-y1)*eta;

        % Gamma2 上的处理
        if (e(5, i)==2)
            VecF=-(tempy-pi*sin(pi*tempy));
            L=zeros(2, 1);
    
            for j1=1:1:2
                temp=VecF.*NE(j1, :)*hxy;
                L(j1)=GaussianQuadrature(0, 1, temp, 4);
            end
    
            % 组装
            for j1=1:1:2
                b(e(j1, i), 1)= b(e(j1, i), 1)+L(j1);
            end
    
        % Gamma3 上的处理
        elseif (e(5, i)==3)
            VecF=-( tempx-pi*sin(pi*tempx));
            L=zeros(2, 1);
        
            for j1=1:1:2
                temp=VecF.*NE(j1, :)*hxy;
                L(j1)=GaussianQuadrature(0, 1, temp, 4);
            end
        
            % 组装
            for j1=1:1:2
                b(e(j1, i), 1)= b(e(j1, i), 1)+L(j1);
            end
        end
    end
    
    % Step 4.4 本质边界条件的处理
    VecD=find(e(5,:)==1 | e(5, :)==4);
    for i=1:1:length(VecD)
        p1=e(1, VecD(i));   p2=e(2, VecD(i));
        A(p1, :)=0; A(p2, :)=0;
        A(p1, p1)=1; A(p2, p2)=1; 
        b(p1)=0; b(p2)=0; 
    end
    
    %% Step 5: 系数矩阵求解
    u=A\b;
    
    %% 画图
    if(ni==NK)
        figure
        pdeplot(model, 'XYData',  u, 'ZData', u);
        figure
        pdemesh(p, e, t, u);
        figure
            pdeplot(p,e,t, "XYData", u);
    end

    %% Step 6: 误差计算
    for i=1:1:length(t)
        x1=p(1, G(1, i));   y1=p(2, G(1, i));
        x2=p(1, G(2, i));   y2=p(2, G(2, i));
        x3=p(1, G(3, i));   y3=p(2, G(3, i));
        S=det([ones(1, 3); [x1, x2, x3]; [y1, y2, y3]])/2;
        DVecy=-[x2-x3, x3-x1, x1-x2];
        DVecx=[y2-y3, y3-y1, y1-y2];
    
        % 求积节点
        tempx=xi*[x1;x2;x3];    tempy=xi*[y1;y2;y3];  
    
        % 真解及其偏导数
        utrue=tempx.*tempy+sin(pi*tempx).*sin(pi*tempy);
        dxu=tempy+pi*cos(pi*tempx).*sin(pi*tempy);
        dyu=tempx+pi*sin(pi*tempx).*cos(pi*tempy);    
    
        % un, dun
        tempu=zeros(Nxi, 1);
        tempdxu=zeros(Nxi, 1);
        tempdyu=zeros(Nxi, 1);
        for j=1:1:3
                tempu=tempu+u(G(j, i))*VecN(:, j);
                tempdxu=tempdxu+u(G(j, i))*VecNx(:, j)*DVecx(j)/(2*S);
                tempdyu=tempdyu+u(G(j, i))*VecNy(:, j)*DVecy(j)/(2*S);
        end
    
        tempu=tempu-utrue;
        tempdxu=tempdxu-dxu;
        tempdyu=tempdyu-dyu;
    
        L2err(ni)=L2err(ni)+S*w*tempu.^2;
        H0err(ni)=H0err(ni)+S*w*(tempdxu.^2+tempdyu.^2);
        H1err(ni)=L2err(ni)+H0err(ni);
        
    end
    
    L2err(ni)=sqrt(L2err(ni));
    H1err(ni)=sqrt(H1err(ni));
    AveH(ni)=sqrt(1/length(t));

end

%% 作图
figure
loglog(AveH, L2err, '--', 'DisplayName', 'L2err');
hold on
loglog(AveH, AveH.^2, '--k', 'DisplayName', 'O(h^2)');
legend
hold off

figure
loglog(AveH, H1err, '--', 'DisplayName', 'H1err');
hold on
loglog(AveH, AveH, '--k', 'DisplayName', 'O(h)');
legend
hold off

%% 数值收敛阶
alpha1=zeros(NK-1, 1);
alpha2=zeros(NK-1, 1);
for i=1:1:NK-1
    alpha1(i)=(log(L2err(i+1))-log(L2err(i)))/(log(AveH(i+1))-log(AveH(i)));
    alpha2(i)=(log(H1err(i+1))-log(H1err(i)))/(log(AveH(i+1))-log(AveH(i)));    
end

