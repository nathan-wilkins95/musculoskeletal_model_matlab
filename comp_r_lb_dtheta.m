function [R_dot, R_ddot, l_bar_dot, l_bar_ddot] = comp_r_lb_dtheta(dR_dtheta, d2R_d2theta, mpartials, theta_dot, theta_ddot)
% function to compute lever dynamics and muscle tendon length dynamics

dl_bar_dtheta       = mpartials(:, 19:20);      
d2l_bar_d2theta     = mpartials(:, 21:22);      
d2lb_dts_dte        = mpartials(:, 23);

% lever dynamics
R_dot  = dR_dtheta*diag(theta_dot);
R_ddot = dR_dtheta*diag(theta_ddot) + d2R_d2theta*diag(theta_dot.^2);

% muscle tendon length dynamics
l_bar_dot   = dl_bar_dtheta*theta_dot;
l_bar_ddot  = d2l_bar_d2theta*(theta_dot.^2) + dl_bar_dtheta*theta_ddot + 2*d2lb_dts_dte*prod(theta_dot); % checked

end