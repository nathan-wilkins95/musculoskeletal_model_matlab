function [ gs, dgs_d_eps_t, d2gs_d2_eps_t ] = comp_gs_pol(eps_t)

global Ws
[ ~,  p ]       = size(Ws);                   p = p - 1;
phi             = zeros(p+1,1);
d_phi_d_theta   = zeros(p+1,1);
d2_phi_d2_theta = zeros(p+1,1);

for i = 1:p+1,      phi(i)             = eps_t^(i-1);               end
for i = 2:p+1,      d_phi_d_theta(i)   = (i-1)*eps_t^(i-2);         end
for i = 3:p+1,      d2_phi_d2_theta(i) = (i-1)*(i-2)*eps_t^(i-3);   end

gs           = Ws*phi;
dgs_d_eps_t   = Ws*d_phi_d_theta;
d2gs_d2_eps_t = Ws*d2_phi_d2_theta;
end