




xmin=0; xmax=1;
N=20; NM=10;

%积分节点
h=(xmax-xmin)/N;
ypoints=xmin+h*(0:1:N);

ErrL=zeros(1, NM);
for i=1:1:NM
    Amatrix=zeros(i, i);
    Bmatrix=zeros(i, i);
    Cmatrix=zeros(N+1, i);
    temp=zeros(i, N+1);
    F=zeros(i, 1);
    
    w=ones(1, N+1)*h;
    w(1)=1.0/2*h; w(N+1)=1.0/2*h;
    for m=1:1:N+1
        for n=1:1:N+1
            temp(m, n)=PolynomialBasis(m, yponits(n));
        end
    end
    Amatrix=temp*diag(w)*temp';
    
    tempBmatrix=zeros(N+1, N+1);
    for m=1:1:N+1
        for n=1:1:N+1
            tempBmatrix(m, n)=-0.5*(ypoints(m)+1)*exp(-ypoints(m)*ypoints(n))*w(n)*w(m);
        end
    end
    Bmatrix=temp*tempBmatrix*temp';

    tempF=zeros(N+1, 1);
    for m=1:1:N+1


    end




end









