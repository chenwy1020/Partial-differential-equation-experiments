%% practice 2: ��1.4.1
%% ȡ���� phi_i=sin(i pi x) �� phi_i=(1-x)x^i ������ 1.4.1

M=50;

err1=zeros(1,M);
err2=zeros(1,M);

condA=zeros(1, M);
condB=zeros(1, M);

x=[0:1:1000]/1000;
u=sin(x)/sin(1)-x;

%% ����1
for i=1:1:M
    N=i;

    %% ��ʼ��
    A=zeros(N, N);
    b=zeros(N, 1);

    %% ���Խ���Ԫ�ز�Ϊ��
    for k=1:1:N
        A(k, k)=((k*pi)^2-1)/2;
        b(k)=-cos(k*pi)/(k*pi);
    end

    %% ���Է��������
    S=A\b;
    
    %% ����L2 ���
    % ���� un ������ڵ㴦ȡֵ 
    C=zeros(1,1001);
    for k=1:1:N
        C=C+S(k)*sin(k*pi*x);
    end
    E=(C-u).^2;

    % ���ø���Simpson��ʽ����L^2 ���
    err1(i)=sqrt(ComSimpson(0, 1, E));
    
    if(i==10) 
       figure
       spy(A);
    end

    %% ������������cond(A, 2) ���� A�Ĵ��� 2 ����
    condA(i)=cond(A, 2);
end



%% ����2
for i=1:1:M
    N=i;
    A=zeros(N, N);
    b=zeros(N, 1);

    %% ����ϵ������
    for k1=1:1:N
        for k2=1:1:N
            A(k1, k2)=k1*k2/(k1+k2-1)-(2*k1*k2+k1+k2)/(k1+k2)+(k1*k2+k1+k2)/(k1+k2+1)...
                +2/(k1+k2+2)-1/(k1+k2+3);
        end
        b(k1)=1/(k1+2)-1/(k1+3);
    end

    %% ���
    T=A\b;
   
    %% ����L2 ���
    % ���� un ������ڵ㴦ȡֵ 
    C=zeros(1,1001);
    for k=1:1:N
        C=C+T(k)*(1-x).*x.^k;
    end
    E=(C-u).^2;

    % ���ø���Simpson��ʽ����L^2 ���
    err2(i)=sqrt(ComSimpson(0, 1, E));

    if(i==10)
       figure
       spy(A);
    end

    %% ������������cond(A, 2) ���� A�Ĵ��� 2 ����
    condB(i)=cond(A, 2 );
end

%% �����ͼ
% loglog() �������ڻ�ָ����������ʾ��仯����
figure
loglog(err1, '-.r', 'DisplayName', 'Basis I');
hold on
loglog(err2, '--k', 'DisplayName', 'Basis II');
legend('Show')
hold off


%% plot the condition number
figure
loglog(condA, '-.r' , 'DisplayName', 'Basis I');
hold on
loglog(condB, '--k', 'DisplayName', 'Basis II');
legend('Show')
hold off