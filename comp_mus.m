function [ ln_dot, Gamma, eps_t, c_alpha, l_bar, mpartials ] = comp_mus( theta, l, q )

global nm L_0 L_t0 Alpha_0

ln_dot = zeros(nm,1);               flg       = zeros(nm,1);
Gamma  = zeros(nm,4);               mpartials = zeros(nm,17);
eps_t  = zeros(nm,1);               c_alpha   = zeros(nm,1);

[ l_bar, dl_bar, d2l_bar ] = comp_l_bar(theta);

mpartials(:,16:17) = [ dl_bar', d2l_bar' ];

for i=1:nm
    % gamma_a, gamma_p
    l_n    = l(i) / L_0(i);
    dl_dln = L_0(i);
    [ gp, ga, dgp, dga ] = comp_g_pa(l_n);
    
    % gamma_s
    c_alpha(i) = sqrt( 1 - ( L_0(i)*sin(Alpha_0(i))/l(i) )^2 );                        % pennation
    eps_t(i)   = ( l_bar(i) - l(i)*c_alpha(i) - L_t0(i)) / L_t0(i);                    % tendon strain
    %
    det_dtheta  = dl_bar(i)/L_t0(i);                                      
    det_dl      = -c_alpha(i)/L_t0(i) - ( ( L_0(i)*sin(Alpha_0(i)))^2 )/( c_alpha(i)*L_t0(i)*l(i)^2 );
    dcos_dl     = ( L_0(i)*sin(Alpha_0(i)) )^2 / ( c_alpha(i)*l(i)^3 );
    d2et_d2l    = -dcos_dl/L_t0(i) + ( ( L_0(i)*sin(Alpha_0(i)))^2 )*( 2*c_alpha(i) + l(i)*dcos_dl )/( L_t0(i)*(l(i)^3)*c_alpha(i)^2 );
    [ gs, dgs, d2gs ] = comp_gs(eps_t(i));
        
    dgs_dtheta = dgs*det_dtheta;
    dgs_dln    = dgs*det_dl*dl_dln;

    % gamma_c
    gc = gs/c_alpha(i) - gp;

    dgc_dtheta = dgs_dtheta / c_alpha(i);
    dgc_dln    = dgs_dln / c_alpha(i) - gs*dcos_dl*dl_dln/c_alpha(i)^2 - dgp*dl_dln;
    
    % gains
    Gamma(i,:) = [ ga, gp, gs, gc ];
    
    [ ln_dot(i), dh, flg(i) ] = ncomp_ln_dot( [ ga, gc ], q(i) );

    mpartials(i,1:15) = [ dgp dga det_dtheta det_dl d2et_d2l dcos_dl dgs d2gs dgs_dtheta dgs_dln dgc_dtheta dgc_dln dh' ];
    
end


