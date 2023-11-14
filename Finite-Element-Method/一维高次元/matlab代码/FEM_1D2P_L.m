%�˳���ΪLagrange��һά��������Ԫ����Ԫ������ͨ�ó���

clear

% ����ѭ���Ĵ���
MK=20;

err2=zeros(MK, 1);      errf=zeros(MK, 1);      
errH1=zeros(MK, 1);   errL2=zeros(MK,1);
CondN=zeros(MK, 1);

for ni=1:1:MK
    n=10*ni;
    a=0; b=1;
    A=zeros(2*n+1, 2*n+1); y=zeros(2*n+1, 1);
    p=zeros(1, 2*n+1); I=zeros(3, n);
    
    %% ������Ϣ
    % ���￼�Ǿ����ʷ� 
    h=(b-a)/(2*n);
    p(1, : )=a+[0:1:2*n]*h; 
    I(1, : )=(1:2:2*n-1);  I(2, : )=(3:2:2*n+1); I(3, : )=(2:2:2*n);

    %% ����������� a(u, v)
    % ���ڴ����� p(x)=1, q(x)=pi^2/4
    p1=1; p2=pi^2/4;

    %% �����׼�����ϵ���ֵ�����ʽ����������ڵ�
    % �������ø���Guass��ʽ
    xi(1)=(a+b)/2-(b-a)/2*sqrt(3/5);
    xi(2)=(a+b)/2;
    xi(3)=(a+b)/2+(b-a)/2*sqrt(3/5);
    
    %% ���㵥Ԫ�նȾ���
    for i=1:1:n
        i1=I(1, i); i2=I(2, i); ic=I(3,i);
        x1=p(i1);  x2=p(i2);  xc=p(ic);
        hi=x2-x1;
    
        % ���� a1, a2, a3, a4
        temp=p1/hi*(4*xi-3).^2+hi*p2*(2*xi-1).^2.*(xi-1).^2;                        a11=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
        temp=4*p1/hi*(4*xi-3).*(1-2*xi)+4*hi*p2*(2*xi-1).*(xi-1).*xi.*(1-xi);       a12=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
        temp=p1/hi*(4*xi-1).*(4*xi-3)+hi*p2*(2*xi-1).^2.*xi.*(xi-1);                a13=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
                                                                                    a21=a12;
        temp=16*p1/hi*(1-2*xi).^2+16*hi*p2*(1-xi).^2.*xi.^2;                        a22=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
        temp=-4*p1/hi*(1-4*xi).*(1-2*xi)+4*hi*p2*(2*xi-1).*xi.^2.*(1-xi);           a23=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
                                                                                    a31=a13;
                                                                                    a32=a23;
        temp=p1/hi*(4*xi-1).^2+hi*p2*(2*xi-1).^2.*xi.^2;                            a33=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0; %���Ժ�a11���
    
        % ���� b1, b2
        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*(1-2*xi).*(-xi);%������ ʵ����Ӧ�����Ҳ࣬��������ϵ��Ϊ�㣬Ҳ�Ͱ���
        b1=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;

        temp=4*hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*xi.*(1-xi);
        bc=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0; 

        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*(2*xi-1).*(xi-1);%�Ҳ����
        b2=(5*temp(1)+8*temp(2)+5*temp(3))/9/2.0;
    
        % ��װ�նȾ���
        A(i1, i1)=A(i1, i1)+a11;  A(i1, ic)=A(i1, ic)+a12;  A(i1, i2)=A(i1, i2)+a13;
        A(ic, i1)=A(i2, i2)+a21;  A(ic, ic)=A(ic, ic)+a22;  A(ic, i2)=A(ic, i2)+a23;
        A(i2, i1)=A(i2, i1)+a31;  A(i2, ic)=A(i2, ic)+a32;  A(i2, i2)=A(i2, i2)+a33;

        y(i1, 1)=y(i1, 1)+b1;     y(ic, 1)=y(ic, 1)+bc;     y(i2, 1)=y(i2, 1)+b2;
    end
    
    % �����ʱ߽�����
    A=A(2:1:2*n+1, 2:1:2*n+1);
    y=y(2:1:2*n+1);

    %% ���ϵ������
    u=zeros(2*n, 1);
    u=A\y;
    u=[0; u];
    
    CondN(ni)=cond(A);

    %% ���� L2 �� H1 ���
    % ���ø���Guass��ʽ������ν�����ֵ����
    % ÿ��������� u=u_{i1} N0(\xi)+ u_{ic} Nc(\xi)  + u_{i2} N1(\xi)
    % ÿ��������� u'=u_{i1} N0'(\xi)+ u_{ic} Nc'(\xi)  + u_{i2} N1'(\xi)
    for i=1:1:n
        i1=I(1, i); i2=I(2, i); ic=I(3,i);
        x1=p(i1);   x2=p(i2);   xc=p(ic);
        hi=x2-x1;

        % u �� ����ֵ
        tempu=u(i1)*(2*xi-1).*(xi-1)+u(ic)*4*xi.*(1-xi)+u(i2)*(1-2*xi).*(-xi);
        tempv=sin(pi/2*(x1+hi*xi));

        % u �� ���ĵ���ֵ
        tempdu=u(i1)*(4*xi-3)/hi+u(ic)*4*(1-2*xi)/hi+u(i2)*(4*xi-1)/hi;
        tempdv=pi/2*cos(pi/2*(x1+hi*xi));

        % ���� L2 ���
        temp=(tempu-tempv).^2*hi/2.0;
        err2(ni)=err2(ni)+(5*temp(1)+8*temp(2)+5*temp(3))/9;

        % ���� H1 �뷶��
        temp=(tempdu-tempdv).^2*hi/2.0;
        errf(ni)=errf(ni)+(5*temp(1)+8*temp(2)+5*temp(3))/9;
    end

    % ���� L2 �� H1 ����
    errH1(ni)=sqrt(errf(ni)+err2(ni));
    errL2(ni)=sqrt(err2(ni));
end

%% ����ֵ��

v=sin(pi/2*p);

figure
plot(p,u,'-.r','DisplayName','u(x)');
title('��ֵ��;�ȷ��ͼ��');
hold on
plot(p,v,'--k','DisplayName','u*(x)');
legend('show');
hold off


%% �����ͼ
% L^2 ���
N=10*[1:1:MK];
h=1./N;
figure
loglog(h, errL2, '--r', 'DisplayName', 'L^2 error');
hold on
loglog(h, h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

% H1 ���
figure
N=10*[1:1:MK];
h=1./N;
loglog(h, errH1, '--r', 'DisplayName', 'H1 error');
hold on
loglog(h, h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

%% Condition number
figure
loglog(h, CondN, '--r', 'DisplayName', 'Condition Number');
hold on
loglog(h, h.^-2, '-.k', 'DisplayName', 'y=h^{-2}');
legend('Show');
hold off


% ��ֵ�׼���
alpha1=zeros(MK-1, 1);
alpha2=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errL2(i+1))-log(errL2(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(errH1(i+1))-log(errH1(i)))/(log(h(i+1))-log(h(i)));    
end

