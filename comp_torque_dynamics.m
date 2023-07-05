function [tau_dot, g_ddot, theta_3dot, M, M_inv, M_dot, M_dot_inv, M_ddot] = comp_torque_dynamics(theta, theta_dot, theta_ddot, tau, R, R_dot, Gamma_s, gamma_s_dot)

global Fmax

%
[M, M_inv]  = comp_M(theta);     [M_dot, M_dot_inv] = comp_M_dot(theta, theta_dot);
M_ddot      = comp_M_ddot(theta, theta_dot, theta_ddot);
[ ~, g ]    = comp_arm_dynamics(theta, theta_dot);
g_dot       = comp_arm_1st_deriv(theta, theta_dot, theta_ddot);
%
tau_dot     = R_dot'*Fmax*Gamma_s + R'*Fmax*gamma_s_dot;
theta_3dot  = M_dot_inv*(tau - g)  + M_inv*(tau_dot - g_dot);
g_ddot      = comp_arm_2nd_deriv(theta, theta_dot, theta_ddot, theta_3dot);
%

end