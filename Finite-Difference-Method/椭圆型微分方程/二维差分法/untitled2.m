clc
clear 
tic

% 进行循环的次数
MK=10;

errC=zeros(MK, 1);    err0=zeros(MK,1);
f  =@(x,y)(2*pi^2*x*y);

for ni=1:1:MK
    n=5*ni; m=4*n^2;
    a=0; b=1;
    u=zeros(n+1, 2*n^2+1);
    px=zeros(n+1,1);    pt=zeros(m+1,1);
    Ix=zeros(n,  2);    It=zeros(m,  2);
    r=1/4;
    
    %% 网格信息
    % 这里考虑均匀剖分 
    hx=(b-a)/n;               ht=(b-a)/m;
    px(:, 1)=a+(0:1:n)*hx;    pt(:, 1)=a+(0:1:m)*ht;
    %Ix(:, 1)=(1:1:n);        Ix(:, 2)=(2:1:n+1); 
    %It(:, 1)=(1:1:m);        It(:, 2)=(2:1:m+1);
    
    u(:,  1)=sin(pi*px);
    %u(1,  :)=0;
    %u(n+1,:)=0;
    for k=1:1:m
        for j=2:1:n
            u(j,k+1)=r*u(j+1,k)+(1-2*r)*u(j,k)+r*u(j-1,k)+f(px(j),pt(k))*ht;
        end
    end


    %% 计算误差
    v=zeros(n+1,1);
    v=exp(-pi^2)*sin(pi*px)+sin(pi*px);

    temp=u(:,2*n^2+1)-v;
    errC(ni)=max(abs(temp));
    err0(ni)=sqrt(sum(temp.^2))*hx;

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
%C 误差
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


