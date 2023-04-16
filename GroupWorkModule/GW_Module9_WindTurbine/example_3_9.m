clear all
close all
clc

% Example 3.9 (Textbook)

a0 = [7.02 0 0]';
a1 = [9.99 4.68 0]';
a2 = [7.44 3.54 2.34]';
a3 = [4.1 2.6 0.6]';

% Input, Output and Feedofrward matrices
B = [0 0 0 1]';
C = [1 0 0 0];
D = 0;

% time, input and initial conditions vector
t = (0:0.01:20)';
u = zeros(size(t,1),1);
x0 = [0 1 1 0]';

for ii = 1:length(a0)
    A = [0 1 0 0; ...
         0 0 1 0; ...
         0 0 0 1; ...
         -a0(ii) -a1(ii) -a2(ii) -a3(ii)];
     
    lambda_A(:,ii) = eig(A);
    
    sys = ss(A,B,C,D);
    [y,t,x] = lsim(sys,u,t,x0);
    
    figure, h1 = subplot(2,1,1); set(h1,'FontName','times','FontSize',16)
    hold on, grid on
    plot(x(:,1),x(:,2),'LineWidth',1)
    xlabel('$$x_1(t)$$','FontName','times','FontSize',16,'Interpreter','latex')
    ylabel('$$x_2(t)$$','FontName','times','FontSize',16,'Interpreter','latex')
    h2 = subplot(2,1,2); set(h2,'FontName','times','FontSize',16)
    hold on, grid on
    plot(t,x(:,1),'k',t,x(:,2),'r',t,x(:,3),'g',t,x(:,4),'m','LineWidth',1)
    xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
    ylabel('$$\mathbf{x}(t)$$','FontName','times','FontSize',16,'Interpreter','latex')
    l1 = legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$');
    set(l1,'FontName','times','FontSize',16,'Interpreter','latex','box','off');
end