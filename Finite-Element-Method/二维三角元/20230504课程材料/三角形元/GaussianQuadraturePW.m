function [ w, xi ] = GaussianQuadraturePW(typeA)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% 对于不同类型的Gauss 公式返回求积节点和权重
% 详见教材 P76 页

if( typeA==1 )
    % xi=zeros(3, 3); w=zeros(1, 3);
    xi=[[2/3, 1/6, 1/6]; [1/6, 2/3, 1/6]; [1/6, 1/6, 2/3]];
    w=1/3*ones(1, 3);
elseif( typeA==2 )
    % xi=zeros(3, 3); w=zeros(1, 3);
    xi=[[0.5, 0.5, 0]; [0.5, 0, 0.5]; [0, 0.5, 0.5]];
    w=1/3*ones(1, 3);
elseif( typeA==3 )
    xi=zeros(4, 3); w=zeros(1, 4);
    xi(1, :)=[1/3, 1/3, 1/3];
    xi(2:4, :)=[[0.6, 0.2, 0.2]; [0.2, 0.6, 0.2]; [0.2, 0.2, 0.6]];
    w=[-0.5625, 0.520833333333333*ones(1,3)];
elseif( typeA==4 )
    xi=zeros(6, 3); w=zeros(1, 6);
    temp=[0.659027622374092, 0.231933368553031, 0.109039009072877];
    xi(1, :)=temp;  xi(2, :)=[temp(1), temp(3), temp(2)];   
    xi(3, :)=[temp(2), temp(1), temp(3)];   xi(4, :)=[temp(2), temp(3), temp(1)];
    xi(5, :)=[temp(3), temp(1), temp(2)];   xi(6, :)=[temp(3), temp(2), temp(1)]; 
    w=1/6*ones(1, 6);
elseif( typeA==5)
    xi=zeros(6, 3); w=zeros(1, 6);
    temp=[0.816847572980459, 0.091576213509771];
    xi(1, :)=[temp(1), temp(2), temp(2)];    xi(2, :)=[temp(2), temp(1), temp(2)];
    xi(3, :)=[temp(2), temp(2), temp(1)];

    temp=[0.108103018168070, 0.445948490915965];
    xi(4, :)=[temp(1), temp(2), temp(2)];    xi(5, :)=[temp(2), temp(1), temp(2)];
    xi(6, :)=[temp(2), temp(2), temp(1)];
    w=[0.109951743655322*ones(1, 3), 0.223381589678011*ones(1, 3)];
elseif( typeA==6 )
    xi=zeros(7, 3); w=zeros(1, 7);
    xi(1, :)=[1/3, 1/3, 1/3];
    temp=[0.736712498968435, 0.237932366472434, 0.025355134551932];
    xi(2, :)=temp;  xi(3, :)=[temp(1), temp(3), temp(2)];   
    xi(4, :)=[temp(2), temp(1), temp(3)];   xi(5, :)=[temp(2), temp(3), temp(1)];
    xi(6, :)=[temp(3), temp(1), temp(2)];   xi(7, :)=[temp(3), temp(2), temp(1)]; 
    w=[0.375, 0.104166666666667*ones(1, 6)];
end

end