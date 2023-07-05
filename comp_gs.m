function [ gamma_s, dgs_det, d2gs_det2 ] = comp_gs(eps_t)
% Compute gamma_s and first and second order derivatives
global eps_th a1 a2 a3

if      eps_t <= 0,         gamma_s   = 0;
                            dgs_det   = 0;
                            d2gs_det2 = 0;
elseif  eps_t <= eps_th,    gamma_s   = a1*( exp(a3*eps_t/eps_th) - 1 ) / ( exp(a3) - 1 );
                            dgs_det   = a1*( exp(a3*eps_t/eps_th)*a3/eps_th  ) / ( exp(a3) - 1 );
                            d2gs_det2 = ( a1/ ( exp(a3) - 1 ))*( exp( (a3/eps_th)*eps_t )  )*(a3/eps_th)^2;
else,                       gamma_s   = a2*( eps_t - eps_th ) + a1;
                            dgs_det   = a2;
                            d2gs_det2 = 0;
            
end

end