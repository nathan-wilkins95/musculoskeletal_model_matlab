function [tau_fe_ddot, T_f, T_e] = comp_tfet(tau, tau_dot, tau_ddot)

global alpha H_f H_e

% c_ddot = -sign(tau).*tau_ddot.*exp( -alpha*abs(tau) )/2 -alpha*sign(tau).*tau_dot.*(-sign(tau).*tau_dot).*exp( -alpha*abs(tau) )/2;
T = [1 1]';

pc_ddot = -sign(tau).*tau_ddot.*exp( -alpha*abs(tau) )/2 -alpha*sign(tau).*tau_dot.*(-sign(tau).*tau_dot).*exp( -alpha*abs(tau) )/2;
qc_ddot = -sign(tau).*T.*exp( -alpha*abs(tau) )/2;% -alpha*sign(tau).*tau_dot.*(-sign(tau).*tau_dot).*exp( -alpha*abs(tau) )/2;

tau_f_ddot  = zeros(2,1);            tau_e_ddot     = zeros(2,1);
T_f         = zeros(2,1);            T_e            = zeros(2,1);            


for i = 1:2
    if tau(i) > 0
        %        
        tau_f_ddot(i)   = tau_ddot(i) + pc_ddot(i);
        tau_e_ddot(i)   = 0 + pc_ddot(i);
        T_f(i)          = T(i) + qc_ddot(i);
        T_e(i)          = 0 + qc_ddot(i);
    else
        %
        tau_f_ddot(i)   = 0 + pc_ddot(i);
        tau_e_ddot(i)   = -tau_ddot(i) + pc_ddot(i);
        T_f(i)          = 0 + qc_ddot(i);
        T_e(i)          = -T(i) + qc_ddot(i);
    end
end

tau_fe_ddot = H_f*tau_f_ddot + H_e*tau_e_ddot;
% Tc          = H_f*T_f + H_e*T_e;                                        

% tau_fe_ddot = [tau_f_ddot(1) -tau_e_ddot(1) tau_f_ddot(2) -tau_e_ddot(2)]';