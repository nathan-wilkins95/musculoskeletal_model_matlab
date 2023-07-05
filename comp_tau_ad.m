function tau = comp_tau_ad(r, q)

global tau_a tau_d

% if r > q
%     tau = tau_a * (0.5 + 1.5*q);
% else
%     tau = tau_d / (0.5 + 1.5*q);
% end

tau = tau_a * (0.5 + 1.5*q);

