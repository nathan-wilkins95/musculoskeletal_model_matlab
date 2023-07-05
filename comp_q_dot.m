function q_dot = comp_q_dot(u, q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute activation dynamics
global tau_a tau_d

%     if u > q,   tau = tau_a * (0.5 + 1.5*q);
%     else,       tau = tau_d / (0.5 + 1.5*q);
%     end
    tau = tau_a * (0.5 + 1.5*q);
    q_dot = (u - q) / tau;
end