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

% Từ đó quy đổi về phương trình rời rạc của đối tượng với T_sample cho
% trước
sys_cont_represent = ss(A_c, B_c, C_c, D_c);
sys_disc_represent = c2d(sys_cont_represent, T_sample, 'zoh');
disp('Kết quả ma trận rời rạc của hệ');
A_d = sys_disc_represent.A
B_d = sys_disc_represent.B
C_d = sys_disc_represent.C
D_d = sys_disc_represent.D

disp('Kích thước biến trạng thái: ')
states_size = size(A_d, 1)

disp('Kích thước ngõ vào: ')
inputs_size = size(B_d, 2)

disp('Kích thước ngõ ra: ')
outputs_size = size(C_d, 1)

Q = eye(states_size);
R = eye(outputs_size);
Q(1,1) = 2;

% Khởi tạo theta
A = A_d
B = B_d
C = C_d
D = D_d
x=[10;2.2];
theta_size = 6;
theta = rand(theta_size,1); % Theta là vector cột
H22 = theta(theta_size);
H21 = [theta(3); theta(5)];
K=inv(H22)*H21';
t=0;
gama=1;
N=10;
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
    H_21 = [theta(3); theta(5)];
    K=H22\(H21');
end
K
[S,K_LQR]=dlqr(A,B,Q,R)