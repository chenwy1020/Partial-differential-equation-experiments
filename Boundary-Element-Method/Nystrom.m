 






xmin=0; xmax=1;
N=(1:1:10)*2;

ErrL=zeros(1, length(N));
hsize=zeros(1,length(N));

for i=1:1:length(N)
    % 权重和求积节点
    Ni=N(i);
    h=(xmax-xmin)/Ni;
    w=ones(1, Ni+1)*h;
    w(1)=1/2*h; w(Ni+1)=1/2*h;
    xpoints=xmin+h*(0:1:Ni);

    % 形成系数矩阵
    Amatrix=zeros(Ni+1, Ni+1);
    F=zeros(Ni+1, 1);
    for m=1:1:Ni+1
        for n=1:1:Ni+1
            Amatrix(m,n)=-0.5*(xpoints(m)+1)*exp(-xpoints(m)*xpoints(n))*w(n);
        end
        F(m)=exp(-xpoints(m))-0.5+0.5*exp(-(xpoints(m)+1));
    end
    M=eye(Ni+1)+Amatrix;

    u=M\F;
    utrue=exp(-xpoints);
    uerr=(utrue'-u).^2;
    ErrL(i)=sqrt(w*uerr);
    hsize(i)=h;
    % 离散L2误差

end

loglog(hsize, ErrL, '-k', 'DisplayName', 'L2 error');
hold on
loglog(hsize,hsize.^2, '--r', 'DisplayName', 'O(h^2)');
legend
hold off
























