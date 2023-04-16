%% TwoTankSystem

clear all
close all
clc

%% Simulation paramters
SIM_TIME = 1000; % simulation time [sec]
STEP_SIZE = 0.05; % fundamental sampling time for running simulations [sec]

%% Model parameters
At = 0.785; % [m^2]
Dv = 2.66; % [m^(1/2)/s]
C0 = 0.056; % [[m^(5/2)/s]]
ka = 0.004; % [m^3/(V*s)]
kh = 2; % [V/m]
kt = 0.1; % [V/C]

% Input and disturbance working point
uw0 = 5; % hot water valve [V]
u10 = uw0;
uc0 = 5; % cold water valve [V]
u20 = uc0;
Av0 = 0.0122; % area of outlet valve [m^2]
v10 = Av0;
Tw0 = 60; % temperature of the hot water [C]
v20 = Tw0;
Tc0 = 30; % temperature of the cold water [C]
v30 = Tc0;

% Simulate the system to evaluate the stationary state in connection with
% the chosen values of input and disturbances
x0 = [1.5,0.2,40,20]'; % initial conditions for the state variables x0 = [H1,H2,T1,T2]'
sim('TwoTankSystem_SteadyState',SIM_TIME,[],[])
x0 = xFinal'; % simulated steady state condition

% Plot of state responses
t = logsout.getElement('H1').Values.Time;
H1 = logsout.getElement('H1').Values.Data;
H2 = logsout.getElement('H2').Values.Data;
Tm1 = logsout.getElement('T1').Values.Data;
Tm2 = logsout.getElement('T2').Values.Data;

figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','Renderer','Painter') 
hf1 = subplot(2,1,1); set(hf1,'FontName','times','FontSize',14)
hold on, grid on
plot(t,H1,'k',t,H2,'r','LineWidth',1.5)
ylabel('[m]','FontName','times','FontSize',14,'Interpreter','latex')
leg1 = legend('$H_1$','$H_2$','Location','NorthEast');
set(leg1,'FontName','times','FontSize',12,'Interpreter','latex','box','off')
hf2 = subplot(2,1,2); set(hf2,'FontName','times','FontSize',14)
hold on, grid on
plot(t,Tm1,'k',t,Tm2,'r','LineWidth',1.5)
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')
ylabel('$[^\circ C]$','FontName','times','FontSize',16,'Interpreter','latex')
leg2 = legend('$T_1$','$T_2$','Location','NorthEast');
set(leg2,'FontName','times','FontSize',12,'Interpreter','latex','box','off')
    
%% Numerically linearize the nonlinear model around the stationary state

% Numerical linearization of a nonlinear system implemented in Simulink can
% be achieved through the command LINMOD (read help)
u0 = [u10 u20 v10 v20 v30]';
[A,B,C,D] = linmod('TwoTankSystem_Linearize',x0,u0);
Bv = B(:,3:5);
B = B(:,1:2);
Dvv = D(:,3:5);
D = D(:,1:2);
x0_lin = zeros(size(A,1),1);

% Present the linearized model
disp('Numerically linearized model of the Two Tank System')
disp('matrix A')
disp(A)
disp('matrix B')
disp(B)
disp('matrix Bv')
disp(Bv)
disp('matrix C')
disp(C)
disp('matrix D')
disp(D)
disp('matrix Dv')
disp(Dvv)

%% Discretize the system

% Choose the sampling time
lambda = eig(A); % eigenvalues
tau = 1./abs(lambda); % time constants

Ts = 1; % sampling time [sec] chosen approximately 10 times smaller than the smallest time constant
[F,G] = c2d(A,B,Ts);
[~,Gv] = c2d(A,Bv,Ts);

% Present the discretized linear model
disp('Discretized linear model of the Two Tank System')
disp('matrix F')
disp(F)
disp('matrix G')
disp(G)
disp('matrix Gv')
disp(Gv)
disp('matrix C')
disp(C)
disp('matrix D')
disp(D)
disp('matrix Dv')
disp(Dv)

% Compare the discretized linear model with the continuous time nonlinear
% model
u0 = [u10 u20 v10 v20 v30]';
v0 = u0(3:5,1);
u0 = u0(1:2,1);
y0 = C*x0;

% Step changes in inputs and disturbances
uw1 = 1.1*uw0;
Tw = 100;
uc1 = 0.9*uc0;
Tc = 400;
Av1 = 1.05*Av0;
Tv = 700;
Tw1 = 1.1*Tw0;
Te1 = 1100;
Tc1 = 1.2*Tc0;
Te2 = 1500;

%% System stability

lambda_OL = eig(A); % open loop eigenvalues
if sum(real(lambda_OL) < 0) == size(A,1)
    disp('Open loop system is asymptotically stable')
elseif sum(real(lambda_OL) < 0) < size(A,1)
    tmp = size(A,1) - sum(real(lambda_OL) < 0);
    if sum(real(lambda_OL) > 0) >= 1 || tmp > 1
        disp('Open loop system is unstable')
    elseif sum(real(lambda_OL) > 0) == 0 && tmp == 1
        disp('Open loop system is marginally stable')
    end
end
TwoTankSys = ss(A,B,C,D);
[w,z] = damp(TwoTankSys);

figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','Renderer','Painter') 
h1 = axes; set(h1,'FontName','times','FontSize',16)
hold on, grid on
h1 = plot(real(lambda_OL),imag(lambda_OL),'x','MarkerSize',15,'MarkerEdgeColor','r','LineWidth',1.5);
line([floor(min(real(lambda_OL))) ceil(max(real(lambda_OL)))+0.1],[0 0],'Color','k','LineWidth',1.5)
line([0 0],[floor(min(imag(lambda_OL)))-1 ceil(max(imag(lambda_OL)))+1],'Color','k','LineWidth',1.5)
xlabel('$\Re(\lambda_i)$','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$\Im(\lambda_i)$','FontName','times','FontSize',16,'Interpreter','latex')
xlim([-0.2 0.1])
ylim([-0.5 0.5])

lambda_CL_des = [-0.1+0.02*1i -0.1-0.02*1i -0.08+0.06*1i -0.08-0.06*1i]'; % desired eigenvalues for the closed loop system
h2 = plot(real(lambda_CL_des),imag(lambda_CL_des),'o','MarkerSize',10,'MarkerEdgeColor','g','MarkerFaceColor','g','LineWidth',1.5);
leg = legend([h1 h2],'Open loop eigenvalues','Desired eigenvalues');
set(leg,'FontName','times','FontSize',14,'Interpreter','latex','box','off')
    
%% Controllability analysis

Mc = ctrb(A,B);
if rank(Mc) == size(A,1)
    disp('Open loop system is controllable')
else
    dimCs = size(A,1) - rank(Mc); % dimension of controllability subspace
    disp('Open loop system is not controllable and the controllability subspace has dimension ',dimCs)
end
    
%% Full state feedback control design (Continuous Time)

K_fs = place(A,B,lambda_CL_des); % continuous time feedback controller gain
K_CT = K_fs;
Ak = A-B*K_fs;
lambda_CL = eig(Ak);
if sum(lambda_CL-lambda_CL_des) < 10^-6
    disp('Closed loop specifications fulfilled')
end

TwoTankSys_CL = ss(Ak,B,C,D);

[y_CL,t_CL] = step(TwoTankSys_CL,150);
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','Renderer','Painter') 
hf1 = subplot(2,2,1); set(hf1,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL,y_CL(:,1,1)/kh,'k','LineWidth',1.5)
ylabel('Water Level $H_2$ [m]','FontName','times','FontSize',14,'Interpreter','latex')
title('$r_1$ [V] Hot Water','FontName','times','FontSize',14,'Interpreter','latex')
hf2 = subplot(2,2,2); set(hf2,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL,y_CL(:,1,2)/kh,'k','LineWidth',1.5)
title('$r_2$ [V] Cold Water','FontName','times','FontSize',14,'Interpreter','latex')
hf3 = subplot(2,2,3); set(hf3,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL,y_CL(:,2,1)/kt,'k','LineWidth',1.5)
ylabel('Temperature $T_2$ $[^\circ C]$','FontName','times','FontSize',16,'Interpreter','latex')
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')
hf4 = subplot(2,2,4); set(hf4,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL,y_CL(:,2,2)/kt,'k','LineWidth',1.5)
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')

%% Set point changes
r0 = [0 0]';
r1 = [0.5 -0.5]';
Tref1 = 100;
Tref2 = 500;

%% Full state feedback + reference feedforward (Continuous Time)

kappa = -C/Ak*B; % steady state gain
N = kappa^-1; % reference feedforward gain
Bff = B*N;

TwoTankSys_CL_FF = ss(Ak,Bff,C,D);

[y_CL_FF,t_CL_FF] = step(TwoTankSys_CL_FF,150);
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','Renderer','Painter') 
hf1 = subplot(2,2,1); set(hf1,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL_FF,y_CL_FF(:,1,1),'k','LineWidth',1.5)
ylabel('Water Level $H_2$ [m]','FontName','times','FontSize',14,'Interpreter','latex')
title('$r_1$ [m] Water Level Tank 2','FontName','times','FontSize',14,'Interpreter','latex')
hf2 = subplot(2,2,2); set(hf2,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL_FF,y_CL_FF(:,1,2),'k','LineWidth',1.5)
title('$r_2$ [$^\circ C$] Temperature Tank 2','FontName','times','FontSize',14,'Interpreter','latex')
hf3 = subplot(2,2,3); set(hf3,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL_FF,y_CL_FF(:,2,1),'k','LineWidth',1.5)
ylabel('Temperature $T_2$ $[^\circ C]$','FontName','times','FontSize',16,'Interpreter','latex')
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')
hf4 = subplot(2,2,4); set(hf4,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL_FF,y_CL_FF(:,2,2),'k','LineWidth',1.5)
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')

%% Set point changes

H_ref0 = y0(1)/kh;
T_ref0 = y0(2)/kt;
H_ref1 = H_ref0 + 0.5;
T_ref1 = T_ref0 + 1;
r0 = [H_ref0 T_ref0]';
Tref1 = 100;
Tref2 = 500;

%% Full state feedback control with integral action (Continuous Time)

lambda_CL_des = [-0.095+0.03*1i -0.095-0.03*1i -0.08+0.06*1i -0.08-0.06*1i -0.05+0.085*1i -0.05-0.085*1i]';
Aof = [A zeros(4,2);-C zeros(2,2)];
Bof = [B;zeros(2,2)];
Br = [zeros(4,2);eye(2)];
Cof = [C zeros(2)];
K_of = place(Aof,Bof,lambda_CL_des)
Ak_of = Aof-Bof*K_of;
K = K_of(:,1:4);
Ki = -K_of(:,5:6);
u0 = u0(1:2,1);
xi0 = inv(Ki)*(u0+K*x0);
xi0_lin = zeros(size(C,1),1);
lambda_CL = eig(Ak_of);
if sum(lambda_CL-lambda_CL_des) < 10^-6
    disp('Closed loop specifications fulfilled')
end
Br(5,1) = kh;
Br(6,2) = kt;
TwoTankSys_CL_i = ss(Ak_of,Br,Cof,zeros(2,2));

[y_CL_i,t_CL_i] = step(TwoTankSys_CL_i);

figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','Renderer','Painter') 
hf5 = subplot(2,2,1); set(hf5,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL_i,y_CL_i(:,1,1)/kh,'k','LineWidth',1.5)
ylabel('Water Level $H_2$ [m]','FontName','times','FontSize',14,'Interpreter','latex')
title('$r_1$ [m] Water Level Tank 2','FontName','times','FontSize',14,'Interpreter','latex')
hf6 = subplot(2,2,2); set(hf6,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL_i,y_CL_i(:,1,2)/kh,'k','LineWidth',1.5)
title('$r_2$ [$^\circ C$] Temperature Tank 2','FontName','times','FontSize',14,'Interpreter','latex')
hf7 = subplot(2,2,3); set(hf7,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL_i,y_CL_i(:,2,1)/kt,'k','LineWidth',1.5)
ylabel('Temperature $T_2$ $[^\circ C]$','FontName','times','FontSize',16,'Interpreter','latex')
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')
hf8 = subplot(2,2,4); set(hf8,'FontName','times','FontSize',14)
hold on, grid on
plot(t_CL_i,y_CL_i(:,2,2)/kt,'k','LineWidth',1.5)
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')

%% Full order observer design

% Observability test
Mo = obsv(A,C);
if rank(Mo) == size(A,1)
    disp('Open loop system is observable')
else
    dimOs = size(A,1) - rank(Mo); % dimension of controllability subspace
    disp('Open loop system is not observable and the observability subspace has dimension ')
    disp(dimOs)
end

lambda_obs_des = 5*lambda_OL;
L = place(A',C',lambda_obs_des)';

lambda_obs = eig(A-L*C);
if sum(lambda_obs-lambda_obs_des) < 10^-6
    disp('Observer dynamics specifications fulfilled')
end

%% Simulate observer

xhat0 = 1.0*x0-x0;
Y0 = C*x0;

% input steps
SIM_TIME = 5000;
time = (0:STEP_SIZE:SIM_TIME)';
Tw = 300; % step time
uw0 = 5;
uw1 = 6.5; % [V]
Tc = 1500; % step time
uc0 = 5;
uc1 = 3.5; % [V]
% disturbance step
Tv = 2200; % step time
Av0 = 0.0122;
Av1 = 1.2*0.0122; % 20% more [m^2]
T1 = 3000; % step time
T2 = 4000; % step time
Tw0 = 60; % [C]
Tw1 = 75; % [C]
Tc0 = 30; % [C]
Tc1 = 20; % [C]

%% Plot 

figure, h1 = axes; set(h1,'FontSize',16,'FontName','times')
hold on, grid on
plot(time,[uw0*ones(length(time(1:Tw/STEP_SIZE)),1); uw1*ones(length(time(Tw/STEP_SIZE+1:end)),1)]./uw0,'LineWidth',1.5)
plot(time,[uc0*ones(length(time(1:Tc/STEP_SIZE)),1); uc1*ones(length(time(Tc/STEP_SIZE+1:end)),1)]./uc0,'LineWidth',1.5)
plot(time,[Av0*ones(length(time(1:Tv/STEP_SIZE)),1); Av1*ones(length(time(Tv/STEP_SIZE+1:end)),1)]./Av0,'LineWidth',1.5)
plot(time,[Tw0*ones(length(time(1:T1/STEP_SIZE)),1); Tw1*ones(length(time(T1/STEP_SIZE+1:end)),1)]./Tw0,'LineWidth',1.5)
plot(time,[Tc0*ones(length(time(1:T2/STEP_SIZE)),1); Tc1*ones(length(time(T2/STEP_SIZE+1:end)),1)]./Tc0,'LineWidth',1.5)
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')
leg1 = legend('$U_w/U_{w0}$','$U_c/U_{c0}$','$A_v/A_{v0}$','$T_w/T_{w0}$','$T_c/T_{c0}$');
set(leg1,'FontName','times','FontSize',14,'Interpreter','latex','Orientation','Horizontal','Location','northoutside')

H1 = logsout{3}.Values.Data(:,1);
H2 = logsout{3}.Values.Data(:,2);
Tp1 = logsout{3}.Values.Data(:,3);
Tp2 = logsout{3}.Values.Data(:,4);
H1_hat = logsout{1}.Values.Data(:,1);
H2_hat = logsout{1}.Values.Data(:,2);
T1_hat = logsout{1}.Values.Data(:,3);
T2_hat = logsout{1}.Values.Data(:,4);

figure, h2 = subplot(2,2,1); set(h2,'FontSize',16,'FontName','times')
hold on, grid on
plot(time,H1,time,H1_hat,'--r','LineWidth',1.5)
ylabel('[m]','FontName','times','FontSize',14,'Interpreter','latex')
leg2 = legend('$H_1$','$\hat{H}_1$');
set(leg2,'FontName','times','FontSize',14,'Interpreter','latex','Orientation','Horizontal')
h3 = subplot(2,2,2); set(h3,'FontSize',16,'FontName','times')
hold on, grid on
plot(time,H2,time,H2_hat,'--r','LineWidth',1.5)
ylabel('[m]','FontName','times','FontSize',14,'Interpreter','latex')
leg3 = legend('$H_2$','$\hat{H}_2$');
set(leg3,'FontName','times','FontSize',14,'Interpreter','latex','Orientation','Horizontal')
h4 = subplot(2,2,3); set(h4,'FontSize',16,'FontName','times')
hold on, grid on
plot(time,Tp1,time,T1_hat,'--r','LineWidth',1.5)
ylabel('[$^\circ C$]','FontName','times','FontSize',14,'Interpreter','latex')
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')
leg4 = legend('$T_1$','$\hat{T}_1$');
set(leg4,'FontName','times','FontSize',14,'Interpreter','latex','Orientation','Horizontal')
h5 = subplot(2,2,4); set(h5,'FontSize',16,'FontName','times')
hold on, grid on
plot(time,Tp2,time,T2_hat,'--r','LineWidth',1.5)
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')
leg5 = legend('$T_2$','$\hat{T}_2$');
set(leg5,'FontName','times','FontSize',14,'Interpreter','latex','Orientation','Horizontal')

%% Changed output equation
C1 = [0 kh 0 0;0 0 kt 0;0 0 0 kt];
Mo1 = obsv(A,C1);
if rank(Mo1) == size(A,1)
    disp('Open loop system is observable')
else
    dimOs = size(A,1) - rank(Mo1); % dimension of controllability subspace
    disp('Open loop system is not observable and the observability subspace has dimension ')
    disp(dimOs)
end

lambda_obs_des = 5*lambda_OL;
L1 = place(A',C1',lambda_obs_des)';

lambda_obs = eig(A-L1*C1);
if sum(lambda_obs-lambda_obs_des) < 10^-6
    disp('Observer dynamics specifications fulfilled')
end

C2 = [0 kh 0 0;0 0 kt 0];
Mo2 = obsv(A,C2);
if rank(Mo2) == size(A,1)
    disp('Open loop system is observable')
else
    dimOs = size(A,1) - rank(Mo2); % dimension of controllability subspace
    disp('Open loop system is not observable and the observability subspace has dimension ')
    disp(dimOs)
end

%% Simulate observer

H1 = logsout{4}.Values.Data(:,1);
H2 = logsout{4}.Values.Data(:,2);
Tp1 = logsout{4}.Values.Data(:,3);
Tp2 = logsout{4}.Values.Data(:,4);
H1_hat = logsout{3}.Values.Data(:,1);
H2_hat = logsout{3}.Values.Data(:,2);
T1_hat = logsout{3}.Values.Data(:,3);
T2_hat = logsout{3}.Values.Data(:,4);

figure, h2 = subplot(2,2,1); set(h2,'FontSize',16,'FontName','times')
hold on, grid on
plot(time,H1,time,H1_hat,'--r','LineWidth',1.5)
ylabel('[m]','FontName','times','FontSize',14,'Interpreter','latex')
leg2 = legend('$H_1$','$\hat{H}_1$');
set(leg2,'FontName','times','FontSize',14,'Interpreter','latex','Orientation','Horizontal')
h3 = subplot(2,2,2); set(h3,'FontSize',16,'FontName','times')
hold on, grid on
plot(time,H2,time,H2_hat,'--r','LineWidth',1.5)
ylabel('[m]','FontName','times','FontSize',14,'Interpreter','latex')
leg3 = legend('$H_2$','$\hat{H}_2$');
set(leg3,'FontName','times','FontSize',14,'Interpreter','latex','Orientation','Horizontal')
h4 = subplot(2,2,3); set(h4,'FontSize',16,'FontName','times')
hold on, grid on
plot(time,Tp1,time,T1_hat,'--r','LineWidth',1.5)
ylabel('[$^\circ C$]','FontName','times','FontSize',14,'Interpreter','latex')
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')
leg4 = legend('$T_1$','$\hat{T}_1$');
set(leg4,'FontName','times','FontSize',14,'Interpreter','latex','Orientation','Horizontal')
h5 = subplot(2,2,4); set(h5,'FontSize',16,'FontName','times')
hold on, grid on
plot(time,Tp2,time,T2_hat,'--r','LineWidth',1.5)
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')
leg5 = legend('$T_2$','$\hat{T}_2$');
set(leg5,'FontName','times','FontSize',14,'Interpreter','latex','Orientation','Horizontal')