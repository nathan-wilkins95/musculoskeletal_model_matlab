function [gamma_s, dgs_det, d2gs_det2] = comp_gs_tpol(eps_t)


gamma_s      = tanh_int(eps_t);
dgs_det      = tanh_fun(eps_t);
d2gs_det2    = tanh_dfun(eps_t);
end