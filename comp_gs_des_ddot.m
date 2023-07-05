function [gs_des_ddot, W] = comp_gs_des_ddot(tau_des, tau_des_dot, R, R_dot, R_ddot)

% gamma_s_des_ddot
tau_fe  = tau_to_tau_fe(tau_des);   tau_fe_dot      = comp_tau_fe_dot(tau_des, tau_des_dot);
W       = comp_w(R);                W_dot           = comp_w_dot(R, R_dot);   
%
W_ddot          = comp_w_ddot(R, R_dot, R_ddot);
% H               = comp_tau_fe_ddot(tau_des, tau_des_dot, tau_des_ddot);    
gs_des_ddot     = 2*W_dot*tau_fe_dot + W_ddot*tau_fe;      % + w*tau_fe_ddot;