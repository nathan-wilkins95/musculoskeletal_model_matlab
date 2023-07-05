function [eps_t, det_dl] = comp_eps_t(t, l, l_bar)

global nm L_0 L_t0 Alpha_0

eps_t = zeros(nm, 1);       c_alpha = zeros(nm, 1);         det_dl = zeros(nm, 1);

for i=1:nm
    c_alpha(i)      = sqrt( 1 - ( L_0(i)*sin(Alpha_0(i))/l(i) )^2 ); 
    eps_t(i)        = ( l_bar(i) - l(i)*c_alpha(i) - L_t0(i)) / L_t0(i); 
    det_dl(i)       = -c_alpha(i)/L_t0(i) - ( ( L_0(i)*sin(Alpha_0(i)))^2 )/( c_alpha(i)*L_t0(i)*l(i)^2 );
end

end