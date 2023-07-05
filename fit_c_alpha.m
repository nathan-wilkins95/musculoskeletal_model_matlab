clear all; close all; clc

params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_alpha shows discontinuity for pennated muscles, can be smoothed

%%
% Define range
ln = [0:0.00001:0.25 0.251:0.000001:0.6 0.601:0.00001:2];

l = L_0*ln;
[~, n] = size(l);

c_alpha = NaN(nm, n);   dcos_dl = NaN(nm, n);
% Compute c_alpha, dcos_dl
for i = 1:n
    [c_alpha(:, i), dcos_dl(:, i)] = comp_alpha(L_0*ln(i));
end

% Define polynomial order for \varPhi
p = 15;
Phi = NaN(p+1,n); d_phi_dl = zeros(p+1,n); 


for i=1:n
    for j = 1:p+1,    Phi(j,i) = ln(i)^(j-1);         end
    for j = 2:p+1,    d_phi_dl(j, i)   = (j-1)*ln(i)^(j-2);         end
end
% Approximate
W = real(c_alpha(1:4, :))*pinv(Phi);

%% Plotting
figure(1)
for i = 1:nm-3
    subplot(nm-3, 1, i), plot(l(i, :), c_alpha(i, :), 'b', l(i, :), W(i, :)*Phi, 'g--'), grid on, ylabel('$\cos (\alpha(l))$ [-]')
end
xlabel('l [m]')

figure(2)
for i = 1:nm-3
    subplot(nm-3, 1, i), plot(l(i, :), dcos_dl(i, :), 'b', l(i, :), W(i, :)*d_phi_dl, 'g--'), grid on, ylabel('$\partial \cos(\alpha) / \partial l$ [-]')
end
xlabel('l [m]')

save('data_files/c_alpha_pol','W');