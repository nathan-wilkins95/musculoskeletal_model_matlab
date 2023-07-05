function X_dot = dynamics_test_phi(t, X)
% dynamics to check projection algorithm (Pomet, Praly)

q = X(1);       q_dot = X(2);

[phi, d_phi_dq] = comp_phi(q);

% projection
if phi <= 0
elseif phi > 0 && d_phi_dq'*q_dot <= 0
else
    q_dot  = q_dot - ( (phi*d_phi_dq*q_dot)/(norm(d_phi_dq)^2) ) * d_phi_dq;
end

q_ddot = -100*q;        % stabilizing term

X_dot = [q_dot q_ddot]';

end