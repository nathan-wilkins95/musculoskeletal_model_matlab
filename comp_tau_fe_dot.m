function tau_fe_dot = comp_tau_fe_dot(tau, tau_dot)

c_dot = comp_c_dot(tau, tau_dot);

tau_f_dot   = zeros(2,1);            tau_e_dot      = zeros(2,1);

for i = 1:2
    if tau(i) > 0
        tau_f_dot(i)  = tau_dot(i) + c_dot(i);
        tau_e_dot(i)  = 0 + c_dot(i);
    else
        tau_f_dot(i)  = 0 + c_dot(i);
        tau_e_dot(i)  = -tau_dot(i) + c_dot(i);
    end
end


tau_fe_dot = [tau_f_dot(1) -tau_e_dot(1) tau_f_dot(2) -tau_e_dot(2)]';