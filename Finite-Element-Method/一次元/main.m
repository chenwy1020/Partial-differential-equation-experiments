clear 
clc

%%
%初始化
N=1000; 

%特别地，取网络分划的中点，避免结点处导数的间断
x=zeros(N,1);
for i=1:1:N
    x(i)=(i-1/2)/N;
end

condA=zeros(20,1);
errL=zeros(20,1);
errH=zeros(20,1);
n=zeros(20,1);
h=zeros(20,1);
du=zeros(N,1);
f=zeros(N,1);
g=zeros(N,1);

for i=1:1:20
    n(i)=10*i;
    h(i)=1/n(i);
end
for i=1:1:20  
    %获取 条件数 数值解的离散点集 和 基函数系数
    [condA(i),u,y]=Finite_element_method(n(i),x);

    %获取 du[N],即 u' 的离散点集
    for j=1:1:N
        k=fix((j*n(i)-1)/N)+1;
        if k==1
            du(j)=y(k)*n(i);
        end
        if k>1
            du(j)=-y(k-1)*n(i)+y(k)*n(i);
        end
    end
    

    %获取 (u*-u)^2 和 (u*-u)^2+(u*'-u')^2 的离散点集
    for j=1:1:N
        f(j)=power(sin(pi*x(j)/2)-u(j),2);
        g(j)=power(cos(pi*x(j)/2)*pi/2-du(j),2)+f(j);
    end

    %获取 L^2误差 和 H^1误差
    errL(i) = sqrt(ComSimpson(0,1,f'));
    errH(i) = sqrt(ComSimpson(0,1,g'));
end

%对比组
v=zeros(N,1);
for i=1:1:N
    v(i)=sin(pi*x(i)/2);
end


figure
plot(x,u,'-.r','DisplayName','u(x)');
title('数值解和精确解图像');
hold on
plot(x,v,'--k','DisplayName','u*(x)');
legend('show');
hold off

% figure
% plot(n,errL,'-.r','DisplayName','errL(n)');
% hold on 
% plot(n,errH,'--k','DisplayName','errH(n)');
% legend('show');
% hold off

figure
plot(log(h),log(errL), '-.r', 'DisplayName', '(log h,log errL)');
legend;

figure
plot(log(h),log(errH), '-.r', 'DisplayName', '(log h,log errH)');
legend;

figure
loglog(h,errL, '-.r', 'DisplayName', 'loglog(errL)');
hold on
loglog(h,errH, '--k', 'DisplayName', 'loglog(errH)');
xlabel('h');
legend('Show');
hold off



figure
loglog(condA,'-.r' , 'DisplayName', 'loglog(condA)');
legend('show');

figure
plot(log(h),log(condA), '-.r', 'DisplayName', '(log h,log errL)');
xlabel('log(h)');
ylabel('log(condA');
legend;




