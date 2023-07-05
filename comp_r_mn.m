function r_mn = comp_r_mn(q, q_dot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute activation dynamics
global tau_a tau_d
    
    if q_dot >= 0,  r_mn = q_dot * tau_a * (0.5 + 1.5*q) + q;
    else,           r_mn = q_dot * tau_d / (0.5 + 1.5*q) + q;
    end
    
end