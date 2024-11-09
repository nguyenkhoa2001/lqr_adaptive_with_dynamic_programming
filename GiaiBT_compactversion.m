clear
clc
% Các tham số của hệ thống
k = 0.05 %Nm
b = 0.1 %N/(m/s)
m = 1 %kg
T_sample = 0.01 % giây

% Khởi tạo ma trận trạng thái của đối tượng
A_c = [0        1;
    -k/m     -b/m];
B_c = [0;
    1/m];
C_c = [1 0];
D_c = 0;

% Từ đó quy đổi về phương trình rời rạc của đối tượng với T_sample
sys_cont_represent = ss(A_c, B_c, C_c, D_c);
sys_disc_represent = c2d(sys_cont_represent, T_sample, 'zoh');
disp('Kết quả ma trận rời rạc của hệ');
A_d = sys_disc_represent.A
B_d = sys_disc_represent.B
C_d = sys_disc_represent.C
D_d = sys_disc_represent.D

states_size = size(A_d, 1);
inputs_size = size(B_d, 2);
outputs_size = size(C_d, 1)

% Khởi tạo theta
A = A_d
B = B_d
C = C_d
D = D_d

x=[0;0.1];

Q = eye(states_size);
Q(1,1) = 2;

R = eye(outputs_size);
theta_size = 6;
theta = rand(theta_size,1); % Theta là vector cột
H22 = theta(theta_size);
H21 = [theta(3); theta(5)];
K=inv(H22)*H21';
t=0;
gama=1;
N=30;
K_data_holder = zeros(2,1,N);
xx=K';
for k=1:N   
    P = eye(theta_size,theta_size);
    Z=[];
    Y=[];
    for i=1:30
        explore=rand;
        xold=x;
        sim('hethong_RR1.mdl');
        phi_1 = [xold(1)^2; xold(1)*xold(2); xold(1)*u; xold(2)^2; xold(2)*u; u^2];
        r=[xold(1);xold(2)]'*Q*[xold(1);xold(2)]+u*R*u;
        u=-K*x;
        phi_2=[x(1)^2;x(1)*x(2);x(1)*u;x(2)^2;x(2)*u;u^2];
        phi=phi_1-phi_2;
        Z=[Z; phi'];
        Y=[Y; r];
    end
    theta=(Z'*Z)\Z'*Y;
    H22=theta(theta_size);
    H21 = [theta(3)/2; theta(5)/2];
    K=H22\(H21');
    K_data_holder(:,:,k) = K;
end
K
[S,K_LQR]=dlqr(A,B,Q,R)

% Extract the data for K1 and K2
K1 = squeeze(K_data_holder(1, 1, :)); % Row 1 data across slices
K2 = squeeze(K_data_holder(2, 1, :)); % Row 2 data across slices
% Define the x-axis (entry indices)
x = 1:N;
% Plot K1 and K2 against the entry index
figure;
plot(x, K1, '-o');hold on;
plot(x, K2, '-s');
hold off;
xlabel('Lần lặp');
ylabel('Giá trị');
title('Đồ thị hội tụ K_1 và K_2', 'Interpreter', 'tex');
legend({'K_1', 'K_2'}, 'Interpreter', 'tex');
grid on;
legend show;
print('K1_K2_plot_entries', '-dsvg');
