function [c_alpha, dcos_dl] = comp_alpha(l)
% Compute pennation angle of a muscle based on its length l
global L_0 Alpha_0 nm
c_alpha = NaN(nm, 1);       dcos_dl = NaN(nm, 1);
    for i = 1:nm
        c_alpha(i) = real(sqrt( 1 - ( L_0(i)*sin(Alpha_0(i))/l(i) )^2 ));
        dcos_dl(i)     = ( L_0(i)*sin(Alpha_0(i)) )^2 / ( c_alpha(i)*l(i)^3 );    % discontinuity candidate
    end
end