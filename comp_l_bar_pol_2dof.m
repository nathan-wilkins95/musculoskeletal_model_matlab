function [ l_bar, dl_bar_dtheta, dl2_bar_d2theta_s, dl2_bar_d2theta_e, dl2_bar_dtheta_s_dtheta_e ] = comp_l_bar_pol_2dof(theta) 

global W6_l_bar

p = 0:1:7; % needs to be the same as in fitting script 
theta_s = theta(1);     theta_e = theta(2);
Phi = comp_vphi(theta_s, theta_e, p)';

npl = length(p);      vphi = zeros(1,npl*npl);

% d_phi_d_theta_s - checked
n = npl + 1;
for i = 2:npl
    a = (p(i))*theta_s^(p(i)-1);
    for j = 1:npl
        b = theta_e^(p(j));
        vphi(n) = a*b;
        n = n+1;
    end
end
d_phi_d_theta_s = vphi';

% d_phi_d_theta_e - checked
vphi = zeros(1, npl*npl); n = 1;
for i = 1:npl
    a = theta_s^(p(i));
    for j = 1:npl
        if j == 1
            vphi(n) = 0;
            n = n+1;
        else
            b = (p(j))*theta_e^(p(j)-1);
            vphi(n) = a*b;
            n = n+1;
        end
    end
end
d_phi_d_theta_e = vphi';

% d2_phi_d2theta_s - checked
vphi = zeros(1,npl*npl);      n = 2*npl+1;
for i = 3:npl
    a = (p(i)-1)*(p(i))*theta_s^(p(i)-2);
    for j = 1:npl
        b = theta_e^(p(j));
        vphi(n) = a*b;
        n = n+1;
    end
end 
d2_phi_d2_theta_s = vphi';

% d2_phi_d2theta_e - checked
vphi = zeros(1,npl*npl);        n = 1;
for i = 1:npl
    a = theta_s^(p(i));
    n = n+2;
    for j = 3:npl
        b = (p(j)-1)*(p(j))*theta_e^(p(j)-2);
        vphi(n) = a*b;
        n = n+1;       
    end
end 
d2_phi_d2_theta_e = vphi';

% mixed term
vphi = zeros(1,npl*npl);        n = npl+1;
for i = 2:npl
    a = (p(i))*theta_s^(p(i)-1);
    n = n+1;
    for j = 2:npl
        b = (p(j))*theta_e^(p(j)-1);
        vphi(n) = a*b;
        n = n+1;
    end
end 
d2_phi_dtheta_s_dtheta_e = vphi';

l_bar                       = W6_l_bar * Phi;
dl_bar_dtheta_s             = W6_l_bar * d_phi_d_theta_s;
dl_bar_dtheta_e             = W6_l_bar * d_phi_d_theta_e;
dl2_bar_d2theta_s           = W6_l_bar * d2_phi_d2_theta_s;
dl2_bar_d2theta_e           = W6_l_bar * d2_phi_d2_theta_e;
dl2_bar_dtheta_s_dtheta_e   = W6_l_bar * d2_phi_dtheta_s_dtheta_e;

dl_bar_dtheta = [dl_bar_dtheta_s dl_bar_dtheta_e];
end