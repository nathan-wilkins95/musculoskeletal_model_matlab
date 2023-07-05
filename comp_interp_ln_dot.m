function [ln_dot, dh] = comp_interp_ln_dot(l_bar, l, q)

global W_ln_dot nm 

nq = 50;        q_r = linspace(0.01, 1, nq);

% normalize l_bar, l
[l_bar_n, l_n] = normalize_l_l_bar(l_bar, l);

pp = 7;      p_vec = 0:1:pp;       
Phi = NaN((pp+1)^2, nm);        ln_dot = NaN(nm, 1);        dh = NaN(nm, 3);

% compute vphi, select W_ln_dot
for m = 1:nm
    Phi(:, m) = comp_vphi(l_bar_n(m), l_n(m), p_vec);  
    
    [qd1, ind1] = min(abs(q_r - q(m)));     % select closed q_r
    if q_r(ind1) <= q(m) % select second q 
        ind2 = ind1 + 1;
    else
        ind2 = ind1;
        ind1 = ind1 - 1;
    end
    q1 = q_r(ind1);         q2 = q_r(ind2);         qd2 = q2 - q(m);
    W1 = W_ln_dot(m, ind1);        W2 = W_ln_dot(m, ind2);
    W1 = W1{1};                    W2 = W2{1};

    Wq = W1  + (W2 - W1) * (q(m) - q1);

    ln_dot(m) = Wq * Phi(:, m);
    
    % compute partials for ln_ddot = W(q)_dot * varPhi + W(q) * varPhi_dot
    dWq_dq = (W2 - W1) * Phi(:, m);
    
    npl = length(p_vec);      vphi = zeros(1,npl*npl);

    % d_phi_d_l_bar - checked
    n = npl + 1;
    for i = 2:npl
        a = (p_vec(i))*l_bar_n(m)^(p_vec(i)-1);
        for j = 1:npl
            b = l_n(m)^(p_vec(j));
            vphi(n) = a*b;
            n = n+1;
        end
    end
    d_phi_d_l_bar = vphi';

    % d_phi_d_l_n - checked
    vphi = zeros(1, npl*npl); n = 1;
    for i = 1:npl
        a = l_bar(m)^(p_vec(i));
        for j = 1:npl
            if j == 1
                vphi(n) = 0;
                n = n+1;
            else
                b = (p_vec(j))*l_n(m)^(p_vec(j)-1);
                vphi(n) = a*b;
                n = n+1;
            end
        end
    end
    d_phi_d_l_n = vphi';
    
    dln_dot_d_l_bar = Wq * d_phi_d_l_bar;
    dln_dot_d_l_n   = Wq * d_phi_d_l_n;
    
    dh(m, :) = [dWq_dq dln_dot_d_l_bar dln_dot_d_l_n];
end

end


