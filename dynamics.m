function [X_dot] = dynamics(t, X)

global nm p8 L_0 Fmax
% Extraction
l           = X(1:nm);      
theta       = X(8:9);     theta_dot   = X(10:11);

% Compute Dynamics
[ ln_dot, Gamma, ~, ~, ~, ~, ~ ] = comp_mus_2dof( theta, l, q );
p8L0    = p8*diag(L_0);     gs      = Gamma(:,3);
%
l_dot  = p8L0*ln_dot;
% tau
[R, ~, ~]   = comp_r_pol(theta); 
Omega       = R'*Fmax;  
tau         = Omega*gs; 
% theta_ddot
[~, M_inv]  = comp_M(theta);
[ ~, gt ]   = comp_arm_dynamics(theta, theta_dot);
theta_ddot  = M_inv*(-gt + tau);
%
X_dot = [l_dot'  theta_dot' theta_ddot']';
end
