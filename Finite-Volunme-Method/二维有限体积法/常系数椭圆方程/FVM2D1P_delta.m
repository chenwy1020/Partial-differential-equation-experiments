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
MK=5;
w=2;

%%声明
errC=zeros(MK, 1);    err0=zeros(MK,1);     err1=zeros(MK,1);   CondN=zeros(MK, 1);
f  =@(x,y)(2*pi^2*x.*y);
dxu=@(y)(y-pi*sin(pi*y));
dyu=@(x)(x-pi*sin(pi*x));
dxv=@(x,y)(y+pi*cos(pi*x)*sin(pi*y));
dyv=@(x,y)(x+pi*sin(pi*x)*cos(pi*y));
h0=0.01*[ w*4, w*3, w*2, w, 1];
h=h0(1:1:MK);

for ni=1:1:MK  

    %% 网格信息
    mesh=generateMesh(model,'GeometricOrder','linear','Hmax',h0(ni));
    [p, e, t ]= meshToPet(mesh);
    C=circumcenter(mesh);
    
    % figure; 
    % pdemesh(p, e, t );

    Np=length(p);
    Nt=length(t);
    Ne=length(e);
    A=zeros(Np,Np);
    b=zeros(Np,1);
    u=zeros(Np,1);
    q=2*pi^2;
    PT=zeros(3, 3);


    %% 组装刚度矩阵
    for i=1:1:Nt
        x1=p(1, t(1 ,i)); y1=p(2, t(1,i));
        x2=p(1, t(2 ,i)); y2=p(2, t(2,i));
        x3=p(1, t(3 ,i)); y3=p(2, t(3,i));
        x0=C(1,i);   y0=C(2,i);
        
        %定义参数矩阵
        PT=triangleTR(x0,y0,x1,y1,x2,y2,x3,y3);

        %形成刚度矩阵
        A(t(1 ,i),t(1 ,i))=A(t(1 ,i),t(1 ,i)) -PT(3,3)/PT(2,1)-PT(3,2)/PT(2,3)+q*PT(1,1);
        A(t(1 ,i),t(2 ,i))=A(t(1 ,i),t(2 ,i)) +PT(3,3)/PT(2,1);
        A(t(1 ,i),t(3, i))=A(t(1 ,i),t(3, i)) +PT(3,2)/PT(2,3);
        b(t(1 ,i))=b(t(1 ,i)) +f(x1,y1)*PT(1,1);

        A(t(2 ,i),t(1 ,i))=A(t(2 ,i),t(1 ,i)) +PT(3,3)/PT(2,1);
        A(t(2 ,i),t(2 ,i))=A(t(2 ,i),t(2 ,i)) -PT(3,3)/PT(2,1)-PT(3,1)/PT(2,2)+q*PT(1,2);
        A(t(2 ,i),t(3, i))=A(t(2 ,i),t(3, i)) +PT(3,1)/PT(2,2);
        b(t(2 ,i))=b(t(2 ,i)) +f(x2,y2)*PT(1,2);
        
        A(t(3 ,i),t(1 ,i))=A(t(3 ,i),t(1 ,i)) +PT(3,2)/PT(2,3);
        A(t(3 ,i),t(2 ,i))=A(t(3 ,i),t(2 ,i)) +PT(3,1)/PT(2,2);
        A(t(3 ,i),t(3, i))=A(t(3 ,i),t(3, i)) -PT(3,2)/PT(2,3)-PT(3,1)/PT(2,2)+q*PT(1,3);
        b(t(3 ,i))=b(t(3 ,i)) +f(x3,y3)*PT(1,3);

    end    
    
    %% 初边值条件处理
    %自然边界条件处理
    for i=1:1:Ne
        x1=p(1, e(1, i));    y1=p(2, e(1, i));
        x2=p(1, e(2, i));    y2=p(2, e(2, i));
        ik1=sqrt((x1-x2)^2+(y1-y2)^2)/2;
        
        if(e(5,i)==2)
            b(e(1,i))=b(e(1,i))-dxu(y1)*ik1;
            b(e(2,i))=b(e(2,i))-dxu(y2)*ik1;
        end

        if(e(5,i)==3)
            b(e(1,i))=b(e(1,i))-dyu(x1)*ik1;
            b(e(2,i))=b(e(2,i))-dyu(x2)*ik1;
        end

    end
    
    
    % 处理本质边界条件
    for i=1:1:Ne
        if(e(5,i)==1)
            A(e(1,i),:)=0;
            A(e(2,i),:)=0;
            A(e(1,i),e(1,i))=1; A(e(2,i),e(2,i))=1;
            b(e(1,i))=0;        b(e(2,i))=0;
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
    for i=1:1:Nt
        x1=p(1, t(1 ,i)); y1=p(2, t(1,i));
        x2=p(1, t(2 ,i)); y2=p(2, t(2,i));
        x3=p(1, t(3 ,i)); y3=p(2, t(3,i));
        S=det([ones(1,3);[x1 x2 x3];[y1 y2 y3]])/2.0;


        tempa=(u(t(1,i))*(y2-y3) +u(t(2,i))*(y3-y1) +u(t(3,i))*(y1-y2))/2/S;
        tempb=(u(t(1,i))*(x3-x2) +u(t(2,i))*(x1-x3) +u(t(3,i))*(x2-x1))/2/S;
        tempf= (tempa-dxv((x1+x2)/2,(y1+y2)/2))^2+(tempb-dyv((x1+x2)/2,(y1+y2)/2))^2 ...
              +(tempa-dxv((x2+x3)/2,(y2+y3)/2))^2+(tempb-dyv((x2+x3)/2,(y2+y3)/2))^2 ...
              +(tempa-dxv((x3+x1)/2,(y3+y1)/2))^2+(tempb-dyv((x3+x1)/2,(y3+y1)/2))^2;
        err1(ni)=err1(ni)+S*tempf/3.0;

    end

    v =(p(1,:).*p(2,:)+sin(pi*p(1,:)).*sin(pi*p(2,:)))';
    temp=u-v;
    errC(ni)=max(abs(temp));
    err0(ni)=sqrt(sum(temp.^2))*h0(ni);
    err1(ni)=sqrt(err1(ni));


end

%% 画数值解
figure
pdeplot(model, 'XYData',  u, 'ZData', u);
figure
pdemesh(p, e, t, u);
figure
pdeplot(p,e,t, "XYData", u);

% %% 画误差图
% C 误差
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
loglog(h, h.^2/4, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

% H1 误差
figure
loglog(h, err1, '--r', 'DisplayName', 'err1');
hold on
loglog(h, h/4, '-.k', 'DisplayName', 'O(h)')
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

function [PT]=triangleTR(x0,y0,x1,y1,x2,y2,x3,y3)
    PT=zeros(3,3);
    S1=det([ones(1,3);[x0 x2 x3];[y0 y2 y3]])/2;
    S2=det([ones(1,3);[x0 x3 x1];[y0 y3 y1]])/2;
    S3=det([ones(1,3);[x0 x1 x2];[y0 y1 y2]])/2;
    PT(1,1)=(S2+S3)/2;
    PT(1,2)=(S3+S1)/2;
    PT(1,3)=(S1+S2)/2;

    PT(2,1)=sqrt((x1-x2)^2+(y1-y2)^2);
    PT(2,2)=sqrt((x2-x3)^2+(y2-y3)^2);
    PT(2,3)=sqrt((x3-x1)^2+(y3-y1)^2);

    PT(3,1)=2*S1/PT(2,2);
    PT(3,2)=2*S2/PT(2,3);
    PT(3,3)=2*S3/PT(2,1);

end

