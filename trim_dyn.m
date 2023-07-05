function meps_dot = trim_dyn(t, meps, theta, tau_des)

global Sigma Fmax nm 

gamma_s = zeros(nm,1);          dgs_det  = zeros(nm,1); 
R = comp_r_pol(theta);

for i = 1:nm, [ gamma_s(i), dgs_det(i), ~ ] = comp_gs_tpol(meps(i)); end
dtau_deps_t = R'*Fmax*diag(dgs_det);

tau = R'*Fmax*gamma_s;

e_tau   = tau_des - tau;
eps_dot = ((dtau_deps_t*Sigma)'*e_tau);


% for i=1:nm
%     if (meps(i) < 0.012) && (eps_dot(i) < 0)                                  % constrained eps_t for DELT1 to not go below 0.0035, the minimum required value for ln_dot to cross 0
%         eps_dot(i) = 0; % constrain with Pomet Praly
%     end
% end
% 
% meps_dot = eps_dot;

% using Pomet/Praly to constrain eps_t
rho_e = 0.027; sigma_e = 0.032;  q_eexp = 4;  eeps_p = 0.55;
eps_t_dot = eps_dot;
for i = 1:nm
    s = abs((meps(i) - rho_e)/sigma_e);

    phi = (2/eeps_p) * ( s.^q_eexp - 1 + eeps_p );

    d_phi_deps = (2/eeps_p) * q_eexp * s.^(q_eexp-2) .* ((meps(i) - rho_e)/sigma_e) * (1/sigma_e);
    
    if phi <= 0
    elseif phi > 0 && d_phi_deps*eps_dot(i) <= 0
    else
        eps_t_dot(i) = - ( (phi*d_phi_deps*eps_dot(i))/(norm(d_phi_deps)^2) ) * d_phi_deps;   
    end
    
end 

meps_dot = eps_t_dot;
end