clear all
clc
close all

% Generation of measured data
SIM_TIME = 1000; % simulation time [s]
STEP_SIZE = 0.001;
Ts = 0.01; % sampling time [s]
time = (0:Ts:SIM_TIME)';
b0 = 0.5*ones(length(time),1); % constant bias [m/s]
v0 = 5*ones(length(time),1); % true velocity [m/s]
v_noise = 0.05*randn(length(time),1); 
v_m = v0 + b0 + v_noise; % measured velocity

x = zeros(length(time),1);
x_noise = randn(length(time),1);

% True position
for ii = 2:length(v0)
    x(ii,1) = x(ii-1,1) + Ts*v0(ii-1,1);
end
x_m = x + x_noise; % measured position

figure, plot(time,v_m)
figure, plot(time,x_m)

% Continuous-time model for the Kalman filter
A=[0 -1;0 0];
B=[1;0];
M=[0;1];
C=[1 0];
D=0;
sysc=ss(A,B,C,D);

% Noise intensities
V=1; % noise intensity of position measurement
q=0.1; W=q^2; % noise intensity of fictitious process noise on bias state

% Compute Discrete-time model
sysd = c2d(sysc,Ts,'zoh');
[F,G,C,D] = ssdata(sysd);

% Noise variances
Vd = V/Ts; % measurement noise variance
Wd_approx = M*W*M'*Ts; % first order approximation from discretization of process noise
Wd=[1/3*q^2*Ts^3 -1/2*q^2*Ts^2;-1/2*q^2*Ts^2 q^2*Ts]; 
% Discrete-time Kalman filter using dlqe:
[Ld1,Pd1,Zd1]=dlqe(F,eye(2),C,Wd,Vd)
[F_kf,G_kf,C_kf,D_kf] = destim(F,G,C,D,Ld1,1,1)
xhat0 = zeros(2,1);
% Discrete-time Kalman filter using lqed (from continuous-time data):
[Ld2,Pd2,Zd2]=lqed(A,M,C,W,V,Ts) %==> same solution !

% Run simulation
sim('Bias_estimation_KF',SIM_TIME);