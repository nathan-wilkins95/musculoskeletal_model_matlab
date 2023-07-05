function [ ln_dot, Gamma, eps_t, c_alpha, l_bar, mpartials, flg ] = comp_mus_2dof( theta, l, q )

global nm L_0 L_t0 Alpha_0

ln_dot = zeros(nm,1);               flg       = zeros(nm,1);
Gamma  = zeros(nm,4);               mpartials = zeros(nm,23);
eps_t  = zeros(nm,1);               c_alpha   = zeros(nm,1);

[ l_bar, dl_bar_dtheta, dl2_bar_d2ts, dl2_bar_d2te, dl2_bar_d_ts_d_te ] = comp_l_bar_pol_2dof(theta);
mpartials(:,19:23) = [ dl_bar_dtheta, dl2_bar_d2ts, dl2_bar_d2te, dl2_bar_d_ts_d_te ];

for i=1:nm
    % gamma_a, gamma_p
    l_n    = l(i) / L_0(i);
    dl_dln = L_0(i);
    [ gp, ga, dgp, dga ] = comp_g_pa(l_n);
    
    % gamma_s
    c_alpha(i)      = sqrt( 1 - ( L_0(i)*sin(Alpha_0(i))/l(i) )^2 );                        % pennation
    eps_t(i)        = ( l_bar(i) - l(i)*c_alpha(i) - L_t0(i)) / L_t0(i);                    % tendon strain
    det_dtheta_s    = dl_bar_dtheta(i, 1) / L_t0(i);        det_dtheta_e    = dl_bar_dtheta(i, 2) / L_t0(i);
    %
    det_dl      = -c_alpha(i)/L_t0(i) - ( ( L_0(i)*sin(Alpha_0(i)))^2 )/( c_alpha(i)*L_t0(i)*l(i)^2 );
    dcos_dl     = ( L_0(i)*sin(Alpha_0(i)) )^2 / ( c_alpha(i)*l(i)^3 );    % discontinuity
    d2et_d2l    = -dcos_dl/L_t0(i) + ( ( L_0(i)*sin(Alpha_0(i)))^2 )*( 2*c_alpha(i) + l(i)*dcos_dl )/( L_t0(i)*(l(i)^3)*c_alpha(i)^2 );

%     [ gs, dgs, d2gs ] = comp_gs(eps_t(i));
%     [ gs, dgs, d2gs ] = comp_gs_pol(eps_t(i));
    [ gs, dgs, d2gs ] = comp_gs_tpol(eps_t(i));
    dgs_dtheta_s = dgs*det_dtheta_s;    dgs_dtheta_e = dgs*det_dtheta_e;
    dgs_dln    = dgs*det_dl*dl_dln;

    % gamma_c
    gc = gs/c_alpha(i) - gp;

    dgc_dtheta_s = dgs_dtheta_s / c_alpha(i);       dgc_dtheta_e = dgs_dtheta_e / c_alpha(i);
    dgc_dln    = dgs_dln / c_alpha(i) - gs*dcos_dl*dl_dln/c_alpha(i)^2 - dgp*dl_dln;
    
    % gains
    Gamma(i,:) = [ ga, gp, gs, gc ];
    % ln_dot
    [ ln_dot(i), dh, flg(i) ] = ncomp_ln_dot( [ ga, gc ], q(i) );

    mpartials(i,1:18) = [ dgp dga det_dtheta_s det_dtheta_e det_dl d2et_d2l dcos_dl dgs d2gs dgs_dtheta_s dgs_dtheta_e dgs_dln dgc_dtheta_s dgc_dtheta_e dgc_dln dh' ];
                        %  1    2       3           4           5       6     7      8    9      10           11         12         13        14         15     16:18    
end


