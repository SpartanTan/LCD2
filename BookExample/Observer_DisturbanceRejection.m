%% Observer design for disturbance rejection
clc
close all
clear all

% System dynamics
zeta = 0.1;
wn1 = 0.3;

A = [0 1;-wn1^2 -2*zeta*wn1];
B = [0 1]';
Bv = [0 0.5]';
C = [1 0];
x0 = [0 0]';

% Disturbance dynamics
disturbance = 'sinusoidal';
switch disturbance
    case 'constant'
        Aw = 0;
        Cw = 1;
        w0 = 0.8;
    case 'sinusoidal'
        wn2 = 0.6;
        Aw = [0 1;-wn2^2 0];
        Cw = [1 0];
        w0 = [0.8 0]';
end

% Augmented system x_a = [x w]'
Axw = Bv*Cw;
Aa = [A Axw;zeros(size(Aw,1),size(A,2)) Aw];
Ba = [B' zeros(1,size(Aw,1))]';
Ca = [C zeros(1,size(Aw,1))];

% Controllability and observability analysis
Mc = ctrb(Aa,Ba);
if rank(Mc) == size(Aa,1)
    disp('System is controllable')
else
    dim_Anc = size(Aa,1)-rank(Mc);
    disp('System uncontrollable subspace has dimension '), disp(dim_Anc) 
    disp('Controllable subspace decomposition')
    [Aac,Bac,Cac,Q,~] = ctrbf(Aa,Ba,Ca)
    Ac = Aac(dim_Anc+1:end,dim_Anc+1:end);
    Bc = Bac(dim_Anc+1:end,1);
    lambda_ctrl_des = [-0.5+0.5*1i,-0.5-0.5*1i];
    Kt1_1 = acker(Ac,Bc,lambda_ctrl_des)
    lambda_ctrl = eig(Ac-Bc*Kt1_1)
    switch disturbance
        case 'constant'
            Kt2_1 = 0; 
            Kt2_2 = 0.5; 
        case 'sinusoidal'
            Kt2_1 = [0 0]; 
            Kt2_2 = [0 0.5]; 
    end
    Kt1 = [Kt2_1 Kt1_1];
    Kt2 = [Kt2_2 Kt1_1];
    K1 = Kt1*Q
    K2 = Kt2*Q
end

Mo = obsv(Aa,Ca);
if rank(Mo) == size(Aa,1)
    disp('System is observable')
    disp('Observer design')
    switch disturbance
        case 'constant'
           lambda_obs_des = [5*lambda_ctrl_des, -2]
           xhat0 = rand(3,1);
        case 'sinusoidal'
           lambda_obs_des = [5*lambda_ctrl_des, -2, -2.35]
           xhat0 = rand(4,1);
    end
    L = acker(Aa',Ca',lambda_obs_des)';
    lambda_obs = eig(Aa-L*Ca)
else
    disp('System observable subspace has dimension '), size(Aa,1)-rank(Mo)
end