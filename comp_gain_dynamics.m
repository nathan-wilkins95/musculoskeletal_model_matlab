function [Eps_t_dot, gamma_s_dot, gamma_p_dot, gamma_a_dot, gamma_c_dot] = comp_gain_dynamics(mpartials, c_alpha, Gamma_s, theta_dot, ln_dot)
% function to compute muscle gain and tendon strain dynamics
global p8 L_0 

dgp_dln             = diag(mpartials(:, 1));    dga_dln             = diag(mpartials(:, 2));    d_eps_t_dtheta      = mpartials(:, 3:4);
d_eps_t_dl          = diag(mpartials(:, 5));    dcos_dl             = mpartials(:, 7);          d_Gamma_s_d_Eps_t   = mpartials(:, 8); 

lm_dot  = p8*diag(L_0)*ln_dot;
%
Eps_t_dot   = d_eps_t_dtheta*theta_dot +  d_eps_t_dl*lm_dot;      % verified
gamma_s_dot = diag(d_Gamma_s_d_Eps_t)*Eps_t_dot;
gamma_p_dot = dgp_dln*ln_dot*p8;
gamma_a_dot = dga_dln*ln_dot*p8;
gamma_c_dot = pinv(diag(c_alpha))*gamma_s_dot - pinv(diag(c_alpha.^2))*diag(Gamma_s)*diag(dcos_dl)*lm_dot - gamma_p_dot; % dcos_dl discontinuous

end