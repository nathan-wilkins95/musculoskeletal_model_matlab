function tau_fe_ddot = comp_tau_fe_ddot(tau, tau_dot, tau_ddot)

c_ddot = comp_c_ddot(tau, tau_dot, tau_ddot);

tau_f_ddot  = zeros(2,1);            tau_e_ddot     = zeros(2,1);

for i = 1:2
    if tau(i) > 0
        %        
        tau_f_ddot(i) = tau_ddot(i) + c_ddot(i);
        tau_e_ddot(i) = 0 + c_ddot(i);
    else
        %
        tau_f_ddot(i) = 0 + c_ddot(i);
        tau_e_ddot(i) = -tau_ddot(i) + c_ddot(i);
    end
end


tau_fe_ddot = [tau_f_ddot(1) -tau_e_ddot(1) tau_f_ddot(2) -tau_e_ddot(2)]';