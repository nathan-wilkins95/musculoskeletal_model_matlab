function [qd_dot] = comp_proj(q_dt, q)
% Proj, equation 12 in Pomet, Praly (1992)
global nm

qd_dot = zeros(nm, 1);

[phi, d_phi_dq] = comp_phi(q);
% qd_dot is 0 if q and q_dot are well behaved
for i = 1:nm
    if phi(i) <= 0
    elseif phi(i) > 0 && d_phi_dq(i)*q_dt(i) <= 0
    else
        qd_dot(i) = - ( (phi(i)*d_phi_dq(i)*q_dt(i))/(norm(d_phi_dq(i))^2) ) * d_phi_dq(i);   
    end
end

end