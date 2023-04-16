%% Example 4.16
% Observer design for marinally stable and unstable system
clear all
clc
close all

%% Marginally stable system
wn = 5; % [rad/s]
A = [0 1;-wn^2 0];
B = [0 1]';
C = [1 0];
x0 = [-1.5 0.25]';

sys_ol = ss(A,B,C,0);
x0 = rand(2,1);

% zero input resposnse
[~,t,x] = initial(sys_ol,x0,10);
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','Renderer','Painter') 
hf1 = axes; set(hf1,'FontName','times','FontSize',14)
hold on, grid on
plot(t,x(:,1),'k',t,x(:,2),'r','LineWidth',1.5)
ylabel('[m]','FontName','times','FontSize',14,'Interpreter','latex')
leg1 = legend('$x_1$','$x_2$','Location','NorthEast');
set(leg1,'FontName','times','FontSize',12,'Interpreter','latex','box','off')
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')

% zero state resposnse
[~,t,x] = step(sys_ol,10);
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','Renderer','Painter') 
hf1 = axes; set(hf1,'FontName','times','FontSize',14)
hold on, grid on
plot(t,x(:,1),'k',t,x(:,2),'r','LineWidth',1.5)
ylabel('[m]','FontName','times','FontSize',14,'Interpreter','latex')
leg1 = legend('$x_1$','$x_2$','Location','NorthEast');
set(leg1,'FontName','times','FontSize',12,'Interpreter','latex','box','off')
xlabel('Time [sec]','FontName','times','FontSize',14,'Interpreter','latex')

lambda_A = eig(A);

%% Oberserver design

wn_obs = 5*wn;
zeta_obs = 0.7;
a = -wn_obs*zeta_obs;
b = sqrt(wn_obs^2-a^2);
lambda_obs_des = [a+b*1i a-b*1i]';

L = acker(A',C',lambda_obs_des)';
lambda_obs = eig(A-L*C);

%% Simulate the observer with the system
xhat0 = x0;
sim('Example_4_16_model')
%%
xhat0 = 5*x0;
sim('Example_4_16_model')

%% Unstable system
clear all

wn_sys = 2;
zeta_sys = -0.3;
a = -zeta_sys*wn_sys;
b = sqrt(wn_sys^2-a^2);

A = [a -b;b a];
B = [0 1]';
C = [1 0];
x0 = [-1.5 0.25]';

lambda_A = eig(A);

wn_obs = 5*wn_sys;
zeta_obs = 0.7;
ao = -wn_obs*zeta_obs;
bo = sqrt(wn_obs^2-a^2);
lambda_obs_des = [ao+bo*1i ao-bo*1i]';

L = acker(A',C',lambda_obs_des)';
lambda_obs = eig(A-L*C);

%% Simulate the observer with the system
xhat0 = x0;
sim('Example_4_16_model')
%%
xhat0 = 5*x0;
sim('Example_4_16_model')