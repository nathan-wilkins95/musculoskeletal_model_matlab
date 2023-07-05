function [ ln_dot ] = comp_mus_ln_dot( l_bar, l, q, l_0, l_t0, alpha_0)     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute ln_dot for single muscle


% gamma_a, gamma_p
l_n    = l / l_0;
[ gp, ga, ~, ~ ] = comp_g_pa(l_n);

% gamma_s
c_alpha      = sqrt( 1 - ( l_0*sin(alpha_0)/l )^2 );                        % pennation
eps_t        = ( l_bar - l*c_alpha - l_t0) / l_t0;                    % tendon strain
%
[ gs, ~, ~ ] = comp_gs_tpol(eps_t);

% gamma_c
gc = gs/c_alpha - gp;

% ln_dot
[ ln_dot, ~, ~ ] = ncomp_ln_dot( [ ga, gc ], q );

end


