%% 测试基函数形状

K1=10; K2=11;
x=[0:1:1000]/1000;

%% 基底一
phi1=sin(K1*pi*x);
phi2=sin(K2*pi*x);

figure
plot(x, phi1, '--r', 'DisplayName', '\phi_1')
hold on
plot(x, phi2, '-.k', 'DisplayName', '\phi_2')
legend('Show')

figure
plot(x, phi1-phi2, '-.k', 'DisplayName', '\phi_1-\phi_2')
legend('Show')


%% 基底二
psi1=(1-x).*x.^K1;
psi2=(1-x).*x.^K2;

dpsi1=K1*x.^(K1-1)-(K1+1)*x.^K1;
dpsi2=K2*x.^(K2-1)-(K2+1)*x.^K2;

figure
plot(x, psi1, '--r', 'DisplayName', '\psi_1')
hold on
plot(x, psi2, '-.k', 'DisplayName', '\psi_2')
legend('Show')

figure
plot(x, psi1-psi2, '-.k', 'DisplayName', '\psi_1-\psi_2')
legend('Show')

plot(x, dpsi1, '--r', 'DisplayName', 'd\psi_1')
hold on
plot(x, dpsi2, '-.k', 'DisplayName', 'd\psi_2')
legend('Show')

figure
plot(x, dpsi1-dpsi2, '-.k', 'DisplayName', 'd\psi_1-d\psi_2')
legend('Show')