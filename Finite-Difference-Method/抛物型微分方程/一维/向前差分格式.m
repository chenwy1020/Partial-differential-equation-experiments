clc
clear 
tic

% 进行循环的次数
MK=10;

err0=zeros(MK,1);
f=@(x,t)(sin(pi*x)+pi^2*t*sin(pi*x));

for ni=1:1:MK
    n=5*ni; m=4*n^2;
    a=0; b=1;
    u=zeros(n+1, m+1);
    px=zeros(n+1,1);    pt=zeros(m+1,1);
    r=1/4;
    
    %% 网格信息
    % 这里考虑均匀剖分 
    hx=(b-a)/n;               ht=(b-a)/m;
    px(:, 1)=a+(0:1:n)*hx;    pt(:, 1)=a+(0:1:m)*ht;
    %Ix(:, 1)=(1:1:n);        Ix(:, 2)=(2:1:n+1); 
    %It(:, 1)=(1:1:m);        It(:, 2)=(2:1:m+1);
    
    u(:,  1)=sin(pi*px);
    %u(:, 1)=0; u(:, m+1)=0;
    for k=1:1:m
        for j=2:1:n
            u(j,k+1)=r*u(j+1,k)+(1-2*r)*u(j,k)+r*u(j-1,k)+f(px(j),pt(k))*ht;
        end
    end


    %% 计算误差
    v=zeros(n+1,1);
    v=exp(-pi^2)*sin(pi*px)+sin(pi*px);

    temp=u(:,m+1)-v;

    err0(ni)=sqrt(temp'*temp*hx);

end

%% 画数值解
V=zeros(n+1,m+1);
[X,T]=meshgrid(0:1/m:1,0:1/n:1);

V=exp(-pi^2*X).*sin(pi*T)+X.*sin(pi*T);

figure
mesh(X,T,u);
hold on;
mesh(X,T,V);
xlabel('x');
ylabel('t');
hold off


%% 画误差图
N=5*(1:1:MK);
h=1./N;
t=1./(4*N.^2);

% H0误差
figure
loglog(t, err0, '--r', 'DisplayName', 'err0');
hold on
loglog(t, t,  '-.k', 'DisplayName', 'O(t)')
legend('Show');
hold off

figure
loglog(h, err0, '--r', 'DisplayName', 'err0');
hold on
loglog(h, 2*h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off
toc
