function [Phi, ln_dot] = comp_vphi_q(nl ,nlb, l_bar, l, q, l_0, l_t0, alpha_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute \varPhi for single q over a grid of l_bar and l
ln_dot = NaN(nlb, nl);

% Normalize
l_bar_n     = (l_bar - min(l_bar)) / ( max(l_bar) - min(l_bar) ); 
l_n_        = (l(:) - min(l(:))) / ( max(l(:)) - min(l(:)) );
l_n         = reshape(l_n_, nlb, nl);

% Define polynomial order
pp = 7;      p_vec = 0:1:pp;      ll = 1;              
Phi = NaN((pp+1)^2, nlb*nl);

for j = 1:nl % muscle length
    for i = 1:nlb % muscle tendon length
        ln_dot(i, j) = comp_mus_ln_dot( l_bar(i), l(i, j), q, l_0, l_t0, alpha_0 );  
        
        Phi(:, ll)  = comp_vphi(l_bar_n(i), l_n(i, j), p_vec); 
        ll = ll + 1;
    end
end