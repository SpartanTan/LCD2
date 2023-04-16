clear all
close all
clc

% Example 3.12 (Textbook)

A = [-3 4 12;...
      1 0 0;...
      0 1 0];

B = [1 0 0]';
C = [1 -1 -2];
D = 0;

% Internal stability
lambda_A = eig(A);

% External stability
[num den] = ss2tf(A,B,C,D);
G = tf(num,den);
[z p k] = tf2zp(num,den);

% Simulate responses
sys = ss(A,B,C,D);
A(1,2) = 4.001;
lambda_Aus = eig(A);
sys_us = ss(A,B,C,D);
[y,t,x] = step(sys);
[y_us,t_us] = step(sys_us);
figure, h1 = subplot(2,1,1); set(h1,'FontName','times','FontSize',16)
hold on, grid on
plot(t,y,'k',t_us,y_us,'r','LineWidth',1)
ax1 = axis;
axis([ax1(1) 5 ax1(3) 0.3]);
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$$y(t)$$','FontName','times','FontSize',16,'Interpreter','latex')
h2 = subplot(2,1,2); set(h2,'FontName','times','FontSize',16)
hold on, grid on
plot(t,x(:,1),'k',t,x(:,2),'r',t,x(:,3),'g','LineWidth',1)
ax2 = axis;
axis([ax2(1) 5 ax2(3) 10]);
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$$x_1(t), \; x_2(t), \; x_3(t)$$','FontName','times','FontSize',16,'Interpreter','latex')