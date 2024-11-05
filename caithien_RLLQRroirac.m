clear
clc

B=[ 0.0450;    0.1038;    0.7378;    0.0007];

A=[0        1.0047      0.0867      -0.0450;  
   -0.0739  0.7490      0.1154      -0.1038;
   -0.5354  -0.3401     0.2303      -0.7378;
   0.0593   0.0316      0.0020      0.9993];

Q=[10 0 0 0;
    0 1 0 0;
    0 0 10 0;
    0 0 0 1];
R=1;
x=[0;0.1;0;0];


xx=[];
uu=[];


%K=[1 1 1 1];
%K=[rand rand rand rand];

% Equal to zeros because we don't need initialize parameters's controller
K=[0 0 0 0];

theta=[rand rand rand rand 2*K(1) rand rand rand 2*K(2) rand rand 2*K(3) rand  2*K(4) 1 ]';
H22=theta(15);
H21=[theta(5)/2; theta(9)/2; theta(12)/2;theta(14)/2];
K=inv(H22).*H21';
K1=K;
t=0;
k=1;
gama=1;
N=4;
xx=K';
tham_so_sai_lech_giua_cac_update = 0.001;
for k=1:N
P=eye(15,15);
    Z=[];
    Y=[];
    Kold = K; % Thầy thêm vào
    for i=1:30
         explore=rand;
         xold=x;
         sim('hethong_RR1.mdl');
         %Phi_1 này là phi lấy bước trước đó
         % và được chọn
         phi_1=[xold(1)^2;xold(1)*xold(2);xold(1)*xold(3);xold(1)*xold(4);xold(1)*u;xold(2)^2;xold(2)*xold(3);xold(2)*xold(4);xold(2)*u;xold(3)^2; xold(3)*xold(4);xold(3)*u;xold(4)^2;xold(4)*u;u^2];
         r=[xold(1);xold(2);xold(3);xold(4)]'*Q*[xold(1);xold(2);xold(3);xold(4)]+u*R*u;
         u=-K*x;
         phi_2=[x(1)^2;x(1)*x(2);x(1)*x(3);x(1)*x(4);x(1)*u;x(2)^2;x(2)*x(3);x(2)*x(4);x(2)*u;x(3)^2;x(3)*x(4); x(3)*u;x(4)^2;x(4)*u;u^2];
         phi=phi_1-phi_2;
     %    m=(phi'*P'*phi);
      %   theta=theta+(P*phi*(r-phi'*theta))/(m);
     %    P=P-(P*phi*phi'*P)/(m);
    %     t=t+1;
          Z=[Z; phi']; 
          Y=[Y; r];
    end
    theta=(Z'*Z)\Z'*Y;
    H22=theta(15);
    H21=[theta(5)/2; theta(9)/2; theta(12)/2;theta(14)/2];
   % S=1/theta(10).*[theta(1) theta(2)/2 theta(3)/2;theta(2)/2 theta(5) theta(6)/2; theta(3)/2 theta(6)/2 theta(8)];
 %   S=[theta(1) theta(2)/2 theta(3)/2;theta(2)/2 theta(5) theta(6)/2; theta(3)/2 theta(6)/2 theta(8)]
    K=1/(H22)*H21';
    if norm(K-Kold) < tham_so_sai_lech_giua_cac_update
        disp("Đã nhỏ hơn sai số đặt ra giữa các lần cập nhật! Thoát")
        break;
    end

end

[S,K_LQR]=dlqr(A,B,Q,R)
%plot(1:N+1,xx(1,:),1:N+1,xx(2,:))