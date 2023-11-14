 






xmin=0; xmax=1;
N=20; NM=10;

%积分节点
h=(xmax-xmin)/N;
ypoints=xmin+h*(0:1:N);

ErrL=zeros(1, NM);
for i=1:1:NM
    CollocationP=(1:1:i)/i*(xmax-xmin);

    Amatrix=zeros(i, i);
    Bmatrix=zeros(i, N+1);
    Cmatrix=zeros(N+1, i);
    F=zeros(i, 1);

    w=ones(1, N+1)*h;
    w(1)=1/2*h; w(N+1)=1/2*h;
    for m=1:1:i
        for n=1:1:i
            Amatrix(m, n)=PolynomialBasis(n, CollocationP(m));
        end
        F(m, 1)=exp(-CollocationP(m))-0.5+0.5*exp(-(CollocationP(m)+1));
    end
   
    for m=1:1:i
        for n=1:1:N+1
            Bmatrix(m, n)=-0.5*(CollocationP(m)+1)*exp(-CollocationP(m)*ypoints(n))*w(n);
            Cmatrix(m, n)=PolynomialBasis(m, ypoints(n));
        end
    end
    
    M=Amatrix+Bmatrix*Cmatrix;
    u=M\F;

    % 离散解
    utrue=exp(-ypoints);
    utest=zeros(1, N+1);
    for m=1:1:i
        utest=utest+u(m)*PolynomialBasis(m, ypoints);
    end
    uerr=(utest-utrue).^2;
    ErrL(i)=sqrt(uerr*w');
end

loglog((1:1:NM), ErrL);



