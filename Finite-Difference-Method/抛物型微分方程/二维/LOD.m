%% LOD
clc
clear 
tic

%% 初始化
m=8;
h=1/m;
t=1/m^2;
a=1/16; r=1;
A1=zeros((m+1)^2,(m+1)^2);
A2=zeros((m+1)^2,(m+1)^2);
b1=zeros((m+1)^2,1);
b2=zeros((m+1)^2,1);
u1=zeros((m+1)^2,1);
u2=zeros((m+1)^2,1);
x=zeros(3,3);

% 初值处理
for i=1:1:m+1
    for j=1:1:m+1
        i1=(i-1)*(m+1)+j;
        u1(i1)=sin(pi*(i-1)*h)*cos(pi*(j-1)*h);
    end
end

% 内点处理
for i=2:1:m
    for j=2:1:m
       i1=(i-1)*(m+1)+j; 
       i2=i1+1;             i3=i1-1;
       i4=i1+m+1;           i5=i1-m-1;
       A1(i1,i1)=1+a*r;     A1(i1,i4)=-a*r/2;    A1(i1,i5)=-a*r/2;
       A2(i1,i1)=1+a*r;     A2(i1,i2)=-a*r/2;    A2(i1,i3)=-a*r/2; 
       
    end
end

% 界点处理
% Neumann边界条件
j=1;
for i=2:1:m
    i1=(i-1)*(m+1)+j;
    i4=i1+m+1;           i5=i1-m-1;   
    A1(i1,i1)=1+a*r;     A1(i1,i4)=-a*r/2;    A1(i1,i5)=-a*r/2;
end
j=m+1;
for i=2:1:m
    i1=(i-1)*(m+1)+j;
    i4=i1+m+1;           i5=i1-m-1;   
    A1(i1,i1)=1+a*r;     A1(i1,i4)=-a*r/2;    A1(i1,i5)=-a*r/2;
end
j=1;
for i=2:1:m
    i1=(i-1)*(m+1)+j;
    i2=i1+1;             
    A2(i1,i1)=1+a*r;     A2(i1,i2)=-a*r;
end
j=m+1;
for i=2:1:m
    i1=(i-1)*(m+1)+j;
    i3=i1-1;             
    A2(i1,i1)=1+a*r;     A2(i1,i3)=-a*r;
end



% Dirichlet边界条件
A1(1:1:m+1,:)=0;   
A2(1:1:m+1,:)=0;  
for i=1:1:m+1
    A1(i,i)=1;
    A2(i,i)=1;
end
A1(m*(m+1)+(1:1:m+1),:)=0; 
A2(m*(m+1)+(1:1:m+1),:)=0; 
for i=1:1:m+1
    A1(m*(m+1)+i,m*(m+1)+i)=1;
    A2(m*(m+1)+i,m*(m+1)+i)=1;
end

%A1=sparse(A1);
%A2=sparse(A2);

%% 按层迭代
for n=1:1:m^2
    % 内点处理(包含Neumann 边界条件)
    for i=2:1:m
        for j=1:1:m+1
            i1=(i-1)*(m+1)+j;
            i4=i1+m+1;           i5=i1-m-1;
            b1(i1)=(1-a*r)*u1(i1)+a*r*u1(i4)/2+a*r*u1(i5)/2;
        end
    end

    % 界点处理
    %Neumann 边界条件
    %Dirichlet 界点条件
    b1(1:1:m+1)=0;
    b1(m*(m+1)+(1:1:m+1))=0;
    
    % 求解 n+1/2层 
    u2=A1\b1;

    % 内点处理(包含Neumann 边界条件)
    for i=2:1:m
        for j=2:1:m
            i1=(i-1)*(m+1)+j;
            i2=i1+1;             i3=i1-1;
            b2(i1)=(1-a*r)*u2(i1)+a*r*u2(i2)/2+a*r*u2(i3)/2;
        end
    end

    % 界点处理
    %Neumann 边界条件
    j=1;
    for i=2:1:m
        i1=(i-1)*(m+1)+j;
        i2=i1+1;
        b2(i1)=(1-a*r)*u2(i1)+a*r*u2(i2);
    end
    j=m+1;
    for i=2:1:m
        i1=(i-1)*(m+1)+j;
        i3=i1-1;
        b2(i1)=(1-a*r)*u2(i1)+a*r*u2(i3);
    end
    % Dirichlet 边界条件
    b2(1:1:m+1)=0;
    b2(m*(m+1)+(1:1:m+1))=0;
    
    % 求解 n+1 层
    u1=A2\b2;

end

n=m/4;
for i=1:1:3
    for j=1:1:3
        i1=n*i+1; i2=n*j+1;
        x(i,j)=u1((i1-1)*(m+1)+i2);
    end
end

toc

