function tau_des = comp_tau_des(Theta, Theta_des)

global control_gains

alpha1 = control_gains(1);              alpha2 = control_gains(2);

theta           = Theta(1:2);        theta_dot          = Theta(3:4);
%
theta_des       = Theta_des(1:2);    theta_des_dot      = Theta_des(3:4); 
theta_des_ddot  = Theta_des(5:6);    

% e1, e2_dot
e1          = theta_des - theta;
theta_s_dot = theta_des_dot + alpha1*e1;
e2          = theta_s_dot - theta_dot;
e1_dot      = - alpha1*e1 + e2;                                            % e1_dot = theta_des_dot - theta_dot;

% theta_ddot
[M, ~]  = comp_M(theta);
[ ~, gt ]   = comp_arm_dynamics(theta, theta_dot);
% tau_des
P           = theta_des_ddot + alpha1*e1_dot + alpha2*e2 + e1;
tau_des     = gt + M*P;
