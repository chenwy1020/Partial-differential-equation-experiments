%% Test by composite Simpson Rule
%% Basis 1: phi_i=sin(i pi x)
%% The integration is computed by compsite Simpson rule with 1001 nodals

M=50;

err1=zeros(1,M);
condA=zeros(1, M);

%% Exact solution
x=[0:1:1000]/1000;
u=sin(x)/sin(1)-x;

for k=1:1:M
    N=k;

    A=zeros(N, N);
    b=zeros(N, 1);

    for i=1:1:N
        tempui=sin(i*x*pi);
        tempdui=i*pi*cos(i*x*pi);
     
        b(i)=ComSimpson(0,1, tempui.*x);
        for j=1:1:N
            tempuj=sin(j*x*pi);
            tempduj=j*pi*cos(j*x*pi);

            A(i, j)=ComSimpson(0,1, tempdui.*tempduj)-ComSimpson(0,1, tempui.*tempuj);
        end
    end

    S=eye(N)/A*b;

    C=zeros(1,1001);
    for i=1:1:N
        C=C+S(i)*sin(i*pi*x);
    end
    E=(C-u).^2;

    % 利用复化Simpson公式计算L^2 误差
    err1(k)=sqrt(ComSimpson(0, 1, E));

    condA(i)=cond(A, 2);
end

loglog([1:1:M], err1, '-.r', 'DisplayName', 'Basis I');
legend('Show')

figure
loglog(condA, '-.r' , 'DisplayName', 'Basis I');
legend('Show');