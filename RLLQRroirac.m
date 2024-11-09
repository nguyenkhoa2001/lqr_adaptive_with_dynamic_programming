clear
clc

B=[ 0.0000;
    0.0100];

A = [1.0000    0.0100;
   -0.0005    0.9990];


Q=[1 0;
    0 1];
R=1;

x=[0;0.1];
T_sample = 0.01;

xx=[];
uu=[];


x_tichluy = [];
theta = rand(6,1);
H22=theta(6);
H21 = [theta(3)/2; theta(5)/2];
K=inv(H22).*H21';
t=0;
gama=1;
N=15;
K_data_holder = zeros(2,1,N);
xx=K';
for k=1:N
P=eye(6,6);
Z=[];
    Y=[];
    for i=1:30
         explore=rand;
         xold=x;
         sim('hethong_RR1.mdl');
         x_tichluy = [x_tichluy, x];
         phi_1 = [xold(1)^2; xold(1)*xold(2); xold(1)*u; xold(2)^2; xold(2)*u; u^2];
         r=[xold(1);xold(2)]'*Q*[xold(1);xold(2)]+u*R*u;
         u=-K*x;
         phi_2=[x(1)^2;x(1)*x(2);x(1)*u;x(2)^2;x(2)*u;u^2];
         phi=phi_1-phi_2;
          Z=[Z; phi']; 
          Y=[Y; r];
    end
    theta=(Z'*Z)\Z'*Y;
    H22=theta(6);
    H21 = [theta(3)/2; theta(5)/2];
    K=H22\(H21')
    K_data_holder(:,:,k) = K;


end

[S,K_LQR]=dlqr(A,B,Q,R)
%plot(1:N+1,xx(1,:),1:N+1,xx(2,:))

K1 = squeeze(K_data_holder(1, 1, :)); % Row 1 data across slices
K2 = squeeze(K_data_holder(2, 1, :)); % Row 2 data across slices
x_plot = 1:N;
figure;
plot(x_plot, K1, '-o');hold on;
plot(x_plot, K2, '-s');
hold off;
xlabel('Lần lặp');
ylabel('Giá trị');
title('Đồ thị hội tụ K_1 và K_2', 'Interpreter', 'tex');
legend({'K_1', 'K_2'}, 'Interpreter', 'tex');
grid on;
legend show;
print('K1_K2_plot_entries', '-dsvg');