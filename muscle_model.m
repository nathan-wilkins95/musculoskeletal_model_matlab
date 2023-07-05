function y = muscle_model(x, l_0 ,alpha_0, l_t0)

global p8
% extraction
l = x(1);       l_bar = x(2);       q = x(3);

% normalized muscle length
l_n = l / l_0;

% passive & active force gains
[ gamma_p, gamma_a ] = comp_g_pa(l_n);

% series force gain
c_alpha = sqrt( 1 - ( l_0*sin(alpha_0)/l )^2 );
eps_t   = ( l_bar - l*c_alpha - l_t0) / l_t0;
gamma_s = comp_gs(eps_t);

% contractile force gain
gamma_c = gamma_s/c_alpha - gamma_p;

% length dynamics
[ ln_dot, mflag ] = comp_ln_dot(gamma_c,gamma_a,q);
lm_dot = p8 * ln_dot * l_0;
% Concatenation
y = [ c_alpha eps_t gamma_p gamma_a gamma_s gamma_c ln_dot mflag lm_dot]';