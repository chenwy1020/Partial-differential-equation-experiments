%% �˳���Ϊһάһ������Ԫ������ͨ�ó���
clear

% ����ѭ���Ĵ���
MK=100;

err2=zeros(MK, 1);      errf=zeros(MK, 1);      
errH1=zeros(MK, 1);   errL2=zeros(MK,1);
CondN=zeros(MK, 1);

for ni=1:1:MK
    n=10+2*ni;
    a=0; b=1;
    A=zeros(n+1, n+1); y=zeros(n+1, 1);
    p=zeros(1, n+1); I=zeros(2, n);
    
    %% ������Ϣ
    % ���￼�Ǿ����ʷ� 
    h=(b-a)/n;
    p=a+[0:1:n]*h;
    I(1, :)=[1:1:n];    I(2, : )=[2:1:n+1];
    
    %% ����������� a(u, v)
    % ���ڴ����� p(x)=1, q(x)=pi^2/4
    p1=1; p2=pi^2/4;

    %% �����׼�����ϵ���ֵ�����ʽ����������ڵ�
    % �������ø���Simpson��ʽ
    M=2; hm=1/(2*M); xi=[0:1:2*M]*hm;
    
    %% ���㵥Ԫ�նȾ���
    for i=1:1:n
        i1=I(1, i); i2=I(2, i);
        x1=p(i1);  x2=p(i2); hi=x2-x1;
    
        % ���� a1, a2, a3, a4
        temp=-p1/hi+hi*p2*(1-xi).*xi;     a1=ComSimpson(0, 1, temp);
        temp=-p1/hi+hi*p2*(1-xi).*xi;     a2=ComSimpson(0, 1, temp);
        temp=p1/hi+hi*p2*xi.^2;   a3=ComSimpson(0, 1, temp);
        temp=p1/hi+hi*p2*(1-xi).^2;   a4=ComSimpson(0, 1, temp);
    
        % ���� b1, b2
        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*(1-xi);
        b1=ComSimpson(0, 1, temp);
    
        temp=hi*pi^2/2*sin(pi/2*(x1+hi*xi)).*xi;
        b2=ComSimpson(0, 1, temp);
    
        % ��װ�նȾ���
        A(i1, i2)=A(i1, i2)+a1;     A(i2, i1)=A(i2, i1)+a2;
        A(i2, i2)=A(i2, i2)+a3;     A(i1, i1)=A(i1, i1)+a4;
        y(i1, 1)=y(i1, 1)+b1;           y(i2, 1)=y(i2, 1)+b2;
    end
    
    % �����ʱ߽���������һ�и�ֵ [1, 0,...]
    A(1, 1)=10^10; 
    y(1)=0;

    %% ���ϵ������
    u=zeros(n+1, 1);
    u=A\y;

    CondN(ni)=cond(A);
    
    %% ���� L2 �� H1 ���
    % ���ø���Simpson��ʽ������ν�����ֵ����
    % ÿ��������� u=u_{i1} N0(\xi)+u_{i2}N1(\xi)
    % ÿ��������� u'=(u_{i2}-u_{i1})/h
    for i=1:1:n
        i1=I(1, i); i2=I(2, i);
        x1=p(i1);  x2=p(i2); hi=x2-x1;

        % u �� ����ֵ
        tempu=u(i1)*(1-xi)+u(i2)*xi;
        tempv=sin(pi/2*(x1+hi*xi));

        % u �� ���ĵ���ֵ
        tempdu=(u(i2)-u(i1))/hi;
        tempdv=pi/2*cos(pi/2*(x1+hi*xi));

        % ���� L2 ���
        temp=(tempu-tempv).^2*hi;
        err2(ni)=err2(ni)+ComSimpson(0, 1, temp);

        % ���� H1 �뷶��
        temp=(tempdu-tempdv).^2*hi;
        errf(ni)=errf(ni)+ComSimpson(0, 1, temp);
    end

    % ���� L2 �� H1 ����
    errH1(ni)=sqrt(errf(ni)+err2(ni));
    errL2(ni)=sqrt(err2(ni));
end

%% �����ͼ
% L^2 ���
N=10+2*[1:1:MK];
h=1./N;
loglog(h, errL2, '--r', 'DisplayName', 'L^2 error');
hold on
loglog(h, h.^2, '-.k', 'DisplayName', 'O(h^2)')
legend('Show');
hold off

% H1 ���
figure
N=10+2*[1:1:MK];
h=1./N;
loglog(h, errH1, '--r', 'DisplayName', 'H1 error');
hold on
loglog(h, h*10, '-.k', 'DisplayName', 'O(h)')
legend('Show');
hold off

%% Condition number
figure
loglog(h, CondN, '--r', 'DisplayName', 'Condition Number');
hold on
loglog(h, h.^-2, '-.k', 'DisplayName', 'y=1/h^2');
legend('Show');
hold off


% ��ֵ�׼���
alpha1=zeros(MK-1, 1);
alpha2=zeros(MK-1, 1);
for i=1:1:MK-1
    alpha1(i)=(log(errL2(i+1))-log(errL2(i)))/(log(h(i+1))-log(h(i)));
    alpha2(i)=(log(errH1(i+1))-log(errH1(i)))/(log(h(i+1))-log(h(i)));    
end

