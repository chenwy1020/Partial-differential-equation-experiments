%% Test by composite Simpson Rule
%% Basis 1: phi_i=(1-x)*x^i
%% The integration is computed by compsite Simpson rule with 1001 nodals

M=50;

err1=zeros(1,M);
condA=zeros(1, M);
errAn=zeros(1, M);
errC=zeros(1, M);

%% Exact solution
x=[0:1:1000]/1000;
u=sin(x)/sin(1)-x;

for k=1:1:M
    N=k;

    A=zeros(N, N);
    b=zeros(N, 1);

    for i=1:1:N
        tempui=(1-x).*x.^i;
        tempdui=i*x.^(i-1)-(i+1)*x.^(i);
     
%         b(i)=ComSimpson(0,1, tempui.*x);
        b(i)=1/(i+2)-1/(i+3);
        for j=1:1:N
            tempuj=(1-x).*x.^j;
            tempduj=j*x.^(j-1)-(j+1)*x.^(j);

            A(i, j)=ComSimpson(0,1, tempdui.*tempduj)-ComSimpson(0,1, tempui.*tempuj);
        end
    end

    EA=zeros(N, N);
    for k1=1:1:N
       for k2=1:1:N
           EA(k1, k2)=k1*k2/(k1+k2-1)-(2*k1*k2+k1+k2)/(k1+k2)+(k1*k2+k1+k2)/(k1+k2+1)...
               +2/(k1+k2+2)-1/(k1+k2+3);
       end
    end

    errAn(k)=norm(EA-A);
    errC(k)=norm(eye(N)/EA*(EA-A)); 

    S=eye(N)/A*b;

    C=zeros(1,1001);
    for i=1:1:N
        C=C+S(i)*(1-x).*x.^i;
    end
    E=(C-u).^2;

    % 利用复化Simpson公式计算L^2 误差
    err1(k)=sqrt(ComSimpson(0, 1, E));

    condA(i)=cond(A, 2);
end

loglog([1:1:M], err1, '-.r', 'DisplayName', 'Basis II');
legend('Show')

figure
loglog(condA, '-.r' , 'DisplayName', 'Basis II');
legend('Show')

figure
plot(x, u, '--r', 'DisplayName', 'u');
hold on
plot(x, C, '--k', 'DisplayName', 'uN');
hold off
legend('Show');