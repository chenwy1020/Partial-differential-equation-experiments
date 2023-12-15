clc
clear 
tic

%% 创建PDE模型
model=createpde;
R1=[3,4,0,1,1,0,0,0,1,1]';
gn=R1;
g=decsg(gn);
geometryFromEdges(model,g);

% 进行循环的次数
MK=6;
a=1; 
%%声明
errC=zeros(MK, 1);    err0=zeros(MK,1);     err1=zeros(MK,1);   CondN=zeros(MK, 1);
f  =@(x,y)(-5*y-2*a);
fv  =@(x,y)(x.*y);
dxu=@(y)(y);
dyu=@(x)(x);
dxv=@(x,y)(y);
dyv=@(x,y)(x);
h0=[0.05, 0.03, 0.02, 0.015, 0.012, 0.01];
%h0=[0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.015, 0.012, 0.01];
h=h0(1:1:MK);
D=@(x,y)([2*x+4,y+a;y+a,x+1]);

for ni=1:1:MK  

    %% 网格信息
    mesh=generateMesh(model,'GeometricOrder','linear','Hmax',h0(ni));
    [p, e, t ]= meshToPet(mesh);
    C=circumcenter(mesh);

    Np=length(p);
    Nt=length(t);
    Ne=length(e);
    A =zeros(Np,Np);
    b =zeros(Np,1);
    u =zeros(Np,1);
    alpha=zeros(2,1,3);
    beta =zeros(2,1,3);
    intD =zeros(2,2,3);
    intf =zeros(3,1);
    S=0;


    %% 组装刚度矩阵
    for i=1:1:Nt
        x1=p(1, t(1 ,i)); y1=p(2, t(1,i));
        x2=p(1, t(2 ,i)); y2=p(2, t(2,i));
        x3=p(1, t(3 ,i)); y3=p(2, t(3,i));
        x0=C(1,i);   y0=C(2,i);
        
        %定义参数矩阵
        [alpha,beta,intD,S,intf]=triangleTR(x0,y0,x1,y1,x2,y2,x3,y3,D,f);

        %形成刚度矩阵
        A(t(1 ,i),t(1 ,i))=A(t(1 ,i),t(1 ,i)) -alpha(:,:,1)'*intD(:,:,1)*beta(:,:,1)/S/2 +alpha(:,:,2)'*intD(:,:,2)*beta(:,:,1)/S/2;
        A(t(1 ,i),t(2 ,i))=A(t(1 ,i),t(2 ,i)) -alpha(:,:,1)'*intD(:,:,1)*beta(:,:,2)/S/2 +alpha(:,:,2)'*intD(:,:,2)*beta(:,:,2)/S/2;
        A(t(1 ,i),t(3, i))=A(t(1 ,i),t(3, i)) -alpha(:,:,1)'*intD(:,:,1)*beta(:,:,3)/S/2 +alpha(:,:,2)'*intD(:,:,2)*beta(:,:,3)/S/2;
        b(t(1 ,i))=b(t(1 ,i)) +intf(1);

        A(t(2 ,i),t(1 ,i))=A(t(2 ,i),t(1 ,i)) -alpha(:,:,3)'*intD(:,:,3)*beta(:,:,1)/S/2 +alpha(:,:,1)'*intD(:,:,1)*beta(:,:,1)/S/2;
        A(t(2 ,i),t(2 ,i))=A(t(2 ,i),t(2 ,i)) -alpha(:,:,3)'*intD(:,:,3)*beta(:,:,2)/S/2 +alpha(:,:,1)'*intD(:,:,1)*beta(:,:,2)/S/2;
        A(t(2 ,i),t(3, i))=A(t(2 ,i),t(3, i)) -alpha(:,:,3)'*intD(:,:,3)*beta(:,:,3)/S/2 +alpha(:,:,1)'*intD(:,:,1)*beta(:,:,3)/S/2;
        b(t(2 ,i))=b(t(2 ,i)) +intf(2);
        
        A(t(3 ,i),t(1 ,i))=A(t(3 ,i),t(1 ,i)) -alpha(:,:,2)'*intD(:,:,2)*beta(:,:,1)/S/2 +alpha(:,:,3)'*intD(:,:,3)*beta(:,:,1)/S/2;
        A(t(3 ,i),t(2 ,i))=A(t(3 ,i),t(2 ,i)) -alpha(:,:,2)'*intD(:,:,2)*beta(:,:,2)/S/2 +alpha(:,:,3)'*intD(:,:,3)*beta(:,:,2)/S/2;
        A(t(3 ,i),t(3, i))=A(t(3 ,i),t(3, i)) -alpha(:,:,2)'*intD(:,:,2)*beta(:,:,3)/S/2 +alpha(:,:,3)'*intD(:,:,3)*beta(:,:,3)/S/2;
        b(t(3 ,i))=b(t(3 ,i)) +intf(3);
    end    
    
    % 初边值条件处理
    % 处理本质边界条件
    for i=1:1:Ne
        x1=p(1, e(1, i));    y1=p(2, e(1, i));
        x2=p(1, e(2, i));    y2=p(2, e(2, i));
        if(e(5,i)==1)
            A(e(1,i),:)=0;
            A(e(2,i),:)=0;
            A(e(1,i),e(1,i))=1; A(e(2,i),e(2,i))=1;
            b(e(1,i))=0;        b(e(2,i))=0;
        end
        if(e(5,i)==2)
            A(e(1,i),:)=0;
            A(e(2,i),:)=0;
            A(e(1,i),e(1,i))=1; A(e(2,i),e(2,i))=1;
            b(e(1,i))=x1*y1;        b(e(2,i))=x2*y2;
        end
        if(e(5,i)==3)
            A(e(1,i),:)=0;
            A(e(2,i),:)=0;
            A(e(1,i),e(1,i))=1; A(e(2,i),e(2,i))=1;
            b(e(1,i))=x1*y1;        b(e(2,i))=x2*y2;
        end
        if(e(5,i)==4)
            A(e(1,i),:)=0;
            A(e(2,i),:)=0;
            A(e(1,i),e(1,i))=1; A(e(2,i),e(2,i))=1;
            b(e(1,i))=0;        b(e(2,i))=0;
        end
    end

    A=sparse(A);

    %% 求解系数矩阵
    u=A\b;

    %% 计算条件数
    CondN(ni)=condest(A);

    %% 计算误差

    % % 计算 L2 和 H1 范数
    v =fv(p(1,:),p(2,:))';

    for i=1:1:Nt
        x1=p(1, t(1 ,i)); y1=p(2, t(1,i));
        x2=p(1, t(2 ,i)); y2=p(2, t(2,i));
        x3=p(1, t(3 ,i)); y3=p(2, t(3,i));
        yi1=(y2+y3)/2; xi1=(x2+x3)/2;
        yj1=(y1+y3)/2; xj1=(x1+x3)/2;
        yk1=(y1+y2)/2; xk1=(x1+x2)/2;

        S=det([ones(1,3);[x1 x2 x3];[y1 y2 y3]])/2.0;

        tempa=(u(t(1,i))*(y2-y3) +u(t(2,i))*(y3-y1) +u(t(3,i))*(y1-y2))/2/S;
        tempb=(u(t(1,i))*(x3-x2) +u(t(2,i))*(x1-x3) +u(t(3,i))*(x2-x1))/2/S;
        tempf= (tempa-dxv((x1+x2)/2,(y1+y2)/2))^2+(tempb-dyv((x1+x2)/2,(y1+y2)/2))^2 ...
              +(tempa-dxv((x2+x3)/2,(y2+y3)/2))^2+(tempb-dyv((x2+x3)/2,(y2+y3)/2))^2 ...
              +(tempa-dxv((x3+x1)/2,(y3+y1)/2))^2+(tempb-dyv((x3+x1)/2,(y3+y1)/2))^2;
        err1(ni)=err1(ni)+S*tempf/3.0;
        tempf= (u(t(1,i))/2+u(t(2,i))/2-fv(xk1,yk1))^2 ...
              +(u(t(2,i))/2+u(t(3,i))/2-fv(xi1,yi1))^2 ...
              +(u(t(3,i))/2+u(t(1,i))/2-fv(xj1,yj1))^2;
        err0(ni)=err0(ni)+S*tempf/3.0;

    end

    temp=u-v;
    errC(ni)=max(abs(temp));
    err0(ni)=sqrt(err0(ni));
    err1(ni)=sqrt(err1(ni));


end

%% 画数值解
figure
pdemesh(p, e, t, u);
hold on
pdemesh(p, e, t, v);
legend('Show');
hold off

% figure
% pdeplot(model, 'XYData',  u, 'ZData', u);
% hold on
% pdeplot(model, 'XYData',  v, 'ZData', v);
% legend('Show');
% hold off

% figure
% pdeplot(p,e,t, "XYData", u);


%% 画误差图
figure
loglog(h, err0, '--* k', 'DisplayName', '||u-u_{h}||_{L^{2}}');
hold on
loglog(h, h.^2/15, '-.g', 'DisplayName', 'O(h^2)');

loglog(h, err1, '--square k', 'DisplayName', '||u-u_{h}||_{H^{1}}');
loglog(h, h/3, '-.b', 'DisplayName', 'O(h)')

loglog(h, errC, '--diamond k', 'DisplayName', '||u-u_{h}||_{L^{\infty}}');
loglog(h, h.^1.8/15, '-.r', 'DisplayName', 'O(h^{1.8})');

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
alpha3=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errC(i+1))-log(errC(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(err0(i+1))-log(err0(i)))/(log(h(i+1))-log(h(i)));    
    alpha3(i)=(log(err1(i+1))-log(err1(i)))/(log(h(i+1))-log(h(i)));
end
toc

function [alpha,beta,intD,S,intf]=triangleTR(x0,y0,x1,y1,x2,y2,x3,y3,D,f)
    w=[5 8 5]'/9.0;
    a=0; b=1;
    I=[-1 0; 0 1];
    t1=(a+b)/2-(b-a)/2*sqrt(3/5.0);
    t2=(a+b)/2;
    t3=(a+b)/2+(b-a)/2*sqrt(3/5.0);
    % 初始化
    alpha=zeros(2,1,3);
    beta =zeros(2,1,3);
    intD =zeros(2,2,3);
    intf=zeros(3,1);

    % 生成各边中点
    yi1=(y2+y3)/2; xi1=(x2+x3)/2;
    yj1=(y1+y3)/2; xj1=(x1+x3)/2;
    yk1=(y1+y2)/2; xk1=(x1+x2)/2;

    % 生成alpha矩阵
    alpha(:,:,1)=I*[yk1-y0; xk1-x0];
    alpha(:,:,2)=I*[yj1-y0; xj1-x0];
    alpha(:,:,3)=I*[yi1-y0; xi1-x0];

    % 生成beta矩阵
    beta(:,:,1)=I*[y3-y2; x3-x2];
    beta(:,:,2)=I*[y1-y3; x1-x3];
    beta(:,:,3)=I*[y2-y1; x2-x1];
    
    % 生成系数矩阵的离散积分公式
    intD(:,:,1)=w(1)*D(x0+(xk1-x0)*t1,y0+(yk1-y0)*t1) ...
               +w(2)*D(x0+(xk1-x0)*t2,y0+(yk1-y0)*t2) ...
               +w(3)*D(x0+(xk1-x0)*t3,y0+(yk1-y0)*t3);
    intD(:,:,2)=w(1)*D(x0+(xj1-x0)*t1,y0+(yj1-y0)*t1) ...
               +w(2)*D(x0+(xj1-x0)*t2,y0+(yj1-y0)*t2) ...
               +w(3)*D(x0+(xj1-x0)*t3,y0+(yj1-y0)*t3);
    intD(:,:,3)=w(1)*D(x0+(xi1-x0)*t1,y0+(yi1-y0)*t1) ...
               +w(2)*D(x0+(xi1-x0)*t2,y0+(yi1-y0)*t2) ...
               +w(3)*D(x0+(xi1-x0)*t3,y0+(yi1-y0)*t3);    
    intD=intD/2;
    
    % 求面积
    S=det([ones(1,3);[x1 x2 x3];[y1 y2 y3]])/2;

    % 求f在Ti的面积分
    S1=det([ones(1,3);[x0 x2 x3];[y0 y2 y3]])/2;
    S2=det([ones(1,3);[x0 x3 x1];[y0 y3 y1]])/2;
    S3=det([ones(1,3);[x0 x1 x2];[y0 y1 y2]])/2;
    tempx0=(x0 +x1 )/2; tempy0=(y0 +y1 )/2;
    tempx1=(xk1+x1 )/2; tempy1=(yk1+y1 )/2;
    tempx2=(x0 +xk1)/2; tempy2=(y0 +yk1)/2;
    tempx3=(x0 +xj1)/2; tempy3=(y0 +yj1)/2;
    tempx4=(xj1+x1 )/2; tempy4=(yj1+y1 )/2;
    intf(1)=(f(tempx0,tempy0)+f(tempx1,tempy1)+f(tempx2,tempy2))*S3/2/3+(f(tempx0,tempy0)+f(tempx3,tempy3)+f(tempx4,tempy4))*S2/2/3;
    tempx0=(x0 +x2 )/2; tempy0=(y0 +y2 )/2;
    tempx1=(xi1+x2 )/2; tempy1=(yi1+y2 )/2;
    tempx2=(x0 +xi1)/2; tempy2=(y0 +yi1)/2;
    tempx3=(x0 +xk1)/2; tempy3=(y0 +yk1)/2;
    tempx4=(xk1+x2 )/2; tempy4=(yk1+y2 )/2;
    intf(2)=(f(tempx0,tempy0)+f(tempx1,tempy1)+f(tempx2,tempy2))*S1/2/3+(f(tempx0,tempy0)+f(tempx3,tempy3)+f(tempx4,tempy4))*S3/2/3;
    tempx0=(x0 +x3 )/2; tempy0=(y0 +y3 )/2;
    tempx1=(xj1+x3 )/2; tempy1=(yj1+y3 )/2;
    tempx2=(x0 +xj1)/2; tempy2=(y0 +yj1)/2;
    tempx3=(x0 +xi1)/2; tempy3=(y0 +yi1)/2;
    tempx4=(xi1+x3 )/2; tempy4=(yi1+y3 )/2;
    intf(3)=(f(tempx0,tempy0)+f(tempx1,tempy1)+f(tempx2,tempy2))*S2/2/3+(f(tempx0,tempy0)+f(tempx3,tempy3)+f(tempx4,tempy4))*S1/2/3;

end

function [intik ,intkj]=triangleTBoundary2(x1,y1,x2,y2,D)
    w=[5 8 5]'/9.0;
    a=0; b=1;
    t1=(a+b)/2-(b-a)/2*sqrt(3/5.0);
    t2=(a+b)/2;
    t3=(a+b)/2+(b-a)/2*sqrt(3/5.0);

    % 生成中点
    yk1=(y1+y2)/2; xk1=(x1+x2)/2;

    f=@(t)([-(yk1-y1), xk1-x1]*D(x1+(xk1-x1)*t,y1+(yk1-y1)*t)*[y1+(yk1-y1)*t, 1]');
    intik=[f(t1),f(t2),f(t3)]*w/2;

    f=@(t)([-(yk1-y2), xk1-x2]*D(x2+(xk1-x2)*t,y2+(yk1-y1)*t)*[y2+(yk1-y2)*t, 1]');
    intkj=-[f(t1),f(t2),f(t3)]*w/2;

end

function [intik ,intkj]=triangleTBoundary3(x1,y1,x2,y2,D)
    w=[5 8 5]'/9.0;
    a=0; b=1;
    t1=(a+b)/2-(b-a)/2*sqrt(3/5.0);
    t2=(a+b)/2;
    t3=(a+b)/2+(b-a)/2*sqrt(3/5.0);

    % 生成中点
    yk1=(y1+y2)/2; xk1=(x1+x2)/2;

    f=@(t)([-(yk1-y1), xk1-x1]*D(x1+(xk1-x1)*t,y1+(yk1-y1)*t)*[1, x1+(xk1-x1)*t]');
    intik=[f(t1),f(t2),f(t3)]*w/2;

    f=@(t)([-(yk1-y2), xk1-x2]*D(x2+(xk1-x2)*t,y2+(yk1-y1)*t)*[1, x2+(xk1-x2)*t]');
    intkj=-[f(t1),f(t2),f(t3)]*w/2;

end

