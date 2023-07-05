function [ r, dr_dtheta, d2r_d2theta ] = comp_r_pol(theta)

global W7_r nm

nj = length(theta);

r = NaN(nm, nj);    dr_dtheta = NaN(nm, nj);    d2r_d2theta = NaN(nm, nj);

for j = 1:length(theta)
    [ ~,  p ]       = size(W7_r{j});                   p = p - 1;
    phi             = zeros(p+1,1);
    d_phi_d_theta   = zeros(p+1,1);
    d2_phi_d2_theta = zeros(p+1,1);

    for i = 1:p+1,      phi(i)             = theta(j)^(i-1);               end
    for i = 2:p+1,      d_phi_d_theta(i)   = (i-1)*theta(j)^(i-2);         end
    for i = 3:p+1,      d2_phi_d2_theta(i) = (i-1)*(i-2)*theta(j)^(i-3);   end

    r(:, j)           = W7_r{j}*phi;
    dr_dtheta(:, j)   = W7_r{j}*d_phi_d_theta;
    d2r_d2theta(:, j) = W7_r{j}*d2_phi_d2_theta;
end
end