% 进行循环的次数
MK=100;

MinEig=zeros(MK, 1);
MaxEig=zeros(MK, 1);
Meshh=zeros(MK, 1);
for ni=1:1:MK
    n=10+2*ni;
    a=0; b=1;
    A=zeros(n+1, n+1); y=zeros(n+1, 1);
    p=zeros(1, n+1); I=zeros(2, n);
    
    %% 网格信息
    % 这里考虑均匀剖分 
    h=(b-a)/n;
    p=a+[0:1:n]*h;
    I(1, :)=[1:1:n];    I(2, : )=[2:1:n+1];
    
    %% 给出变分问题 a(u, v)
    % 对于此问题 p(x)=1, q(x)=pi^2/4
    p1=1; p2=-2.0001;

    %% 定义标准区域上的数值求积公式，定义求积节点
    % 这里利用复化Simpson公式
    M=10; hm=1/(2*M); xi=[0:1:2*M]*hm;
    
    %% 计算单元刚度矩阵
    for i=1:1:n
        i1=I(1, i); i2=I(2, i);
        x1=p(i1);  x2=p(i2); hi=x2-x1;
    
        % 计算 a1, a2, a3, a4
        temp=-p1/hi+hi*p2*(1-xi).*xi;     a1=ComSimpson(0, 1, temp);
        temp=-p1/hi+hi*p2*(1-xi).*xi;     a2=ComSimpson(0, 1, temp);
        temp=p1/hi+hi*p2*xi.^2;   a3=ComSimpson(0, 1, temp);
        temp=p1/hi+hi*p2*(1-xi).^2;   a4=ComSimpson(0, 1, temp);
    
        % 计算 b1, b2
        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*(1-xi);
        b1=ComSimpson(0, 1, temp);
    
        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*xi;
        b2=ComSimpson(0, 1, temp);
    
        % 组装刚度矩阵
        A(i1, i2)=A(i1, i2)+a1;     A(i2, i1)=A(i2, i1)+a2;
        A(i2, i2)=A(i2, i2)+a3;     A(i1, i1)=A(i1, i1)+a4;
        y(i1, 1)=y(i1, 1)+b1;           y(i2, 1)=y(i2, 1)+b2;
    end

    [s, t]=eig(A);
    t=diag(t);
    MinEig(ni)=min(abs(t));
    MaxEig(ni)=max(abs(t));
    Meshh(ni)=h;
end

loglog(Meshh, MaxEig, '-r', 'DisplayName', 'MaxEig');
hold on
loglog(Meshh, Meshh.^-1, '--k', 'DisplayName', 'h');
legend('Show')

figure
loglog(Meshh, MinEig, '-r', 'DisplayName', 'MinEig');
hold on
loglog(Meshh, Meshh, '--k', 'DisplayName', 'h');
legend('Show')
