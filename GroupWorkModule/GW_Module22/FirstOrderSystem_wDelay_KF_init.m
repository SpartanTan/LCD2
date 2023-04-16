%% Kalman filtering of first order system with time delay

clc
close all
clear al

SIM_TIME = 20;
STEP_SIZE = 0.01;

% System parameter
alpha = -1; % [rad/s]
A = alpha;
B = 1;
Bw = 1;
C = 1;
D = 0;
x0 = 20;

tau = [0 0.1 0.5 1]; % time delay [s]
T = tau(4);

% Measurement noise
Vn = 1; % measurement noise intensity
Tc = 0.01; % correlation time of the noise

% Kalman filter design 
A_kf = A;
B_kf = B;
Bw_kf = Bw;
C_kf = C;
D_kf = D;

Vw = 100; % process noise intensity

[L,P1,lambda] = lqe(A_kf,Bw_kf,C_kf,Vw,Vn);
xhat0 = 0;

% Generate measurements
sim('FirstOrderSystem_wDelay_KF',SIM_TIME);

time = logsout.getElement(1).Values.Time;
n = logsout.getElement(1).Values.Data;
x = logsout.getElement(3).Values.Data;
xhat = logsout.getElement(2).Values.Data;
u_tau = logsout.getElement(5).Values.Data;
u = logsout.getElement(6).Values.Data;

y = x+n*sqrt(Tc);

figure, h1 = axes; set(h1,'FontName','times','FontSize',16)
hold on, grid on
plot(time,x,'k','LineWidth',2)
plot(time,xhat,'r','LineWidth',2)
plot(time,y,'g','LineWidth',1)
plot(time,u,'m','LineWidth',0.5)
plot(time,u_tau,'--c','LineWidth',0.5)

%% Include 2nd order Padé approximation in the model of the Kalman filter
[numP,denP] = pade(T,2);
[Atau,Btau,Ctau,Dtau] = tf2ss(numP,denP);

A_kf = [A B*Ctau;zeros(2,1) Atau];
B_kf = [B*Dtau; Btau];
Bw_kf = [Bw; zeros(2,1)];
C_kf = [C zeros(1,2)];
D_kf = 0;

Vw = 1;

[L,P2,lambda] = lqe(A_kf,Bw_kf,C_kf,Vw,Vn);
xhat0 = zeros(3,1);

% Generate measurements
sim('FirstOrderSystem_wDelay_KF',SIM_TIME);

time = logsout.getElement(1).Values.Time;
n = logsout.getElement(1).Values.Data;
x = logsout.getElement(3).Values.Data;
xhat = logsout.getElement(2).Values.Data;
u_tau = logsout.getElement(5).Values.Data;
u = logsout.getElement(6).Values.Data;

y = x+n*sqrt(Tc);

figure, h1 = axes; set(h1,'FontName','times','FontSize',16)
hold on, grid on
plot(time,x,'k','LineWidth',2)
plot(time,xhat,'r','LineWidth',2)
plot(time,y,'g','LineWidth',1)
plot(time,u,'m','LineWidth',0.5)
plot(time,u_tau,'--c','LineWidth',0.5)
