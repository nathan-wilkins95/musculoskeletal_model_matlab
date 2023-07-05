function meps_dot = trim_dyn_1dof(t, meps, theta_e)

global Sigma Fmax nm W7_e

gamma_s = zeros(nm,1);          dgs_det  = zeros(nm,1); 
r_e       = comp_r_pol(theta_e, W7_e);

for i = 1:nm, [ gamma_s(i), dgs_det(i) ] = comp_gs(meps(i)); end

% dtau_s_deps_t = r_s'*Fmax*diag(dgs_det);
% dtau_e_deps_t = r_e'*Fmax*diag(dgs_det);

dtau_deps_t = r_e'*Fmax*diag(dgs_det);

tau_d = 0;
tau_e   = r_e'*Fmax*gamma_s;

e_tau = tau_d - tau_e;
meps_dot = ((dtau_deps_t*Sigma)'*e_tau);