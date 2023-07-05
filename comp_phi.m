function [phi, d_phi_dq] = comp_phi(q)
% equation 10 in Pomet, Praly (1992)
global eps_p q_exp rho sigma

s = abs((q - rho)/sigma);

phi = (2/eps_p) * ( s.^q_exp - 1 + eps_p );

d_phi_dq = (2/eps_p) * q_exp * s.^(q_exp-2) .* ((q - rho)/sigma) * (1/sigma);

end



