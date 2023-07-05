close all, clear, clc

params

%% Define Space

nt = 150;        nl = 100;        nq = 50;  % sampling frequencies
eps_t_max = 0.05;      eps_t_min = -0.001; % maximum value for eps_t

% define ranges 
eps_t_r = linspace(-0.001, eps_t_max, nl);      
eps_t   = zeros(nt, nl);                        q       = linspace(0.01, 1, nq);

for i = 1:nl, eps_t(:, i) = eps_t_r(i)*ones(nt, 1);   end

% l_bar
l_bar_struct    = load('data_files/l_bar_data.mat');    % import l_bar values from OpenSim
muscle_data     = load('data_files/m_specs.mat');
m_specs = muscle_data;      muscle_names = fieldnames(m_specs);     nlb = 100;

% theta_s     = linspace(-pi/2, pi/2, nlb);
% theta_e     = linspace(-deg2rad(5), 2.7, nlb);
L_bar         = zeros(nm, nlb*nlb);     L_bar_      = zeros(nm, nt);

for m = 1:nm 
    muscle_name = muscle_names(m);
    data = l_bar_struct.(string(muscle_name));
    L_bar(m, :) = reshape(data.', 1, []); % make sure to flatten by rows
    L_bar_(m, :) = linspace(min(L_bar(m, :)), max(L_bar(m, :)), nt);       % turn l_bar into linear space
end

l = NaN(nt, nl, nm);        gamma_s = NaN(nt, nl);  
n =  101;                   tw_init = linspace(0,1e-3,n);       tic
% compute corresponding l
for j = 1:nl
    for i = 1:nt 
        lt = L_t0*eps_t(i, j) + L_t0;
%         l(i, j, :) = (L_bar_(:, i) - L_t0*eps_t(i, j) - L_t0); % TODO: needs to be adjusted for pennated muscles, 
        % namely probably define l and then cut off depending on eps_t 
        l0 = L_0;
        
        [ t_i2, l_i ] = ode45(@(t,y) trim_length(t, y, L_bar_(:, i), lt), tw_init, l0, options);
        l(i, j, :) = l_i(end, :);

%         gamma_s(i, j) = comp_gs(eps_t(i, j));  % needs to be changed for
%         more muscles
    end
end
toc


% 
% get min max values for l
% l_min = NaN(nm, 1);     l_max = NaN(nm, 1);
% for m = 1:nm
%     l_m = squeeze(l(:, :, m));
%     l_min(m) = min(l_m(:));
%     l_max(m) = max(l_m(:));
% end
% l_minmax = [l_min l_max]
% save('data_files/l_minmax', "l_minmax")
%% Plot eps_t, gamma_s
figure(1)
for i = 1:nm
        subplot(2, 4, i), surf(repmat(L_bar_(i, :)', 1, nl), squeeze(l(:, :, i)), eps_t, 'EdgeAlpha',0), grid on
%         subplot(2, 4, i), surf(repmat(L_bar_(i, :)', 1, nl), squeeze(l(:, :, i)), squeeze(Eps_t(:, :, i)), 'EdgeAlpha',0), grid on
    xlabel('$\bar l$ [deg]'), ylabel('$l$ [m]'), zlabel('$\epsilon_{\rm t}$ [-]')
end


% figure(11)
% surf(repmat(L_bar_(1, :)', 1, nl), l, gamma_s, 'EdgeAlpha',0), grid on
% xlabel('$\theta$ [deg]'), ylabel('$l$ [m]'), zlabel('$\gamma_{\rm s}$ [-]')

% figure(12)
% surf(repmat(L_bar_(1, :)', 1, nl), eps_t, gamma_s, 'EdgeAlpha',0), grid on
% xlabel('$\bar l$ [deg]'), ylabel('$\epsilon_{\rm t}$ [-]'), zlabel('$\gamma_{\rm s}$ [-]')
%% Generate Data

PHI = cell(nm, nq);      LN_DOT = cell(nm, nq);

for i = 1:nq
    for m = 1:nm
        [PHI{m, i}, LN_DOT{m, i}] = comp_vphi_q(nl, nt, L_bar_(m, :), squeeze(l(:, :, m)), q(i), L_0(m), L_t0(m), Alpha_0(m));
    end
end


%% Approximation 
      
W = cell(nm, nq);        LN_DOT_HAT = cell(nm, nq);
for m = 1:nm
    disp(muscle_names{m})
    for i = 1:nq
    
    
        ln_dot = LN_DOT{m, i};        Phi = PHI{m, i};
        ln_dot_v = ln_dot(:);         % ln_dot_vn = (ln_dot_v - min(ln_dot_v)) ./ ( max(ln_dot_v) - min(ln_dot_v));


        W_ln_dot = ln_dot_v' * pinv(Phi);   % ((Phi' * Phi + 1e-09*eye(n*nl))\Phi'); %

        ln_dot_v_hat = W_ln_dot*Phi;       % Ln_dot_hat = reshape(ln_dot_v_hat', nt, nl);     

        disp("\hat{\dot l_{\rm n}} - \dot l_{\rm n} = "+num2str(mean((ln_dot_v_hat - ln_dot_v').^2)/max(ln_dot_v)))
%         W_ln_dot
        W{m, i} = W_ln_dot;        
        LN_DOT_HAT{m, i} = ln_dot_v_hat;
    end

end
% save model
save("data_files/ln_dot_pol", 'W')

%% Plotting
nth = 10;

% Theta = repmat(theta', 1, nl);      theta_v = Theta(:);     
eps_t_v = eps_t(:);

for m = 1:nm
figure(m+1)       
l_bar_v = L_bar_(m, :);        L_bar_r = repmat(l_bar_v', 1, nl);       l_bar_rv = L_bar_r(:);
for i = 1:nq
ln_dot_v_hat = LN_DOT_HAT{m, i};   Ln_dot = LN_DOT{m, i};
subplot(4, 5, i),   surf(repmat(l_bar_v', 1, nl), eps_t, Ln_dot, 'EdgeAlpha', 0, 'FaceAlpha',0.8), grid on, hold on
                    scatter3(l_bar_rv(1:nth:end), eps_t_v(1:nth:end), ln_dot_v_hat(1:nth:end), 'red','filled'), hold off
end
xlabel('$\bar l$ [m]'), ylabel('$\epsilon_{\rm t}$ [-]'), zlabel('$\dot l_{\rm n}$ [1/s]')
end

% i = 1; m = 4;
% ln_dot_v_hat = LN_DOT_HAT{m, i};   Ln_dot = LN_DOT{m, i};
% figure(100)
% surf(repmat(l_bar_v', 1, nl), eps_t, Ln_dot, 'EdgeAlpha', 0, 'FaceAlpha',0.8), grid on, hold on
%                     scatter3(l_bar_rv(1:nth:end), eps_t_v(1:nth:end), ln_dot_v_hat(1:nth:end), 'red','filled'), hold off
%                     xlabel('$\bar l$ [m]'), ylabel('$\epsilon_{\rm t}$ [-]'), zlabel('$\dot l$ [m/s]')
%%
% figure(2)
% surf(repmat(theta', 1, nl)/d2r, l, squeeze(gamma_s_dot(:, :, 2))), grid on
% xlabel('$\theta$ [deg]'), ylabel('$l$ [m]'), zlabel('$\dot \gamma_{\rm s}$ [-]')

% figure(22)
% surf(repmat(theta', 1, nl)/d2r, eps_t, squeeze(gamma_s_dot(:, :, 2))), grid on
% xlabel('$\theta$ [deg]'), ylabel('$\epsilon_{\rm t}$ [-]'), zlabel('$\dot \gamma_{\rm s}$ [-]')

% figure(23)
% for i = 1:nq
% ln_dot_v_hat = LN_DOT_HAT{i};   Ln_dot = LN_DOT{i};
% subplot(4, 5, i),   surf(repmat(theta', 1, nl)/d2r, eps_t, Ln_dot, 'EdgeAlpha', 0, 'FaceAlpha',0.8), grid on, hold on
%                     scatter3(theta_v(1:nth:end)/d2r, eps_t_v(1:nth:end), ln_dot_v_hat(1:nth:end), 'red','filled'), hold off
% end
% xlabel('$\theta$ [deg]'), ylabel('$\epsilon_{\rm t}$ [-]'), zlabel('$\dot l$ [m/s]')

% figure(24)
% surf(repmat(theta', 1, nl)/d2r, l, Ln_dot, 'EdgeAlpha', 0), grid on
% xlabel('$\theta$ [deg]'), ylabel('$\l $ [m]'), zlabel('$\dot l$ [m/s]')
% 
% figure(25)
% surf(repmat(theta', 1, nl)/d2r, l, Ln_dot_hat, 'EdgeAlpha', 0), grid on
% xlabel('$\theta$ [deg]'), ylabel('$\l $ [m]'), zlabel('$\hat{\dot l}$ [1/s]')

% figure(26)
% surfc(repmat(theta', 1, nl)/d2r, eps_t, Ln_dot_hat, 'EdgeAlpha', 0), grid on
% xlabel('$\theta$ [deg]'), ylabel('$\epsilon_{\rm t} $ [m]'), zlabel('$\hat{\dot l}$ [1/s]')
% 
% figure(27)
% surfc(repmat(theta', 1, nl)/d2r, eps_t, Ln_dot_hat-Ln_dot, 'EdgeAlpha',0), grid on
% xlabel('$\theta$ [deg]'), ylabel('$\epsilon_{\rm t} $ [m]'), zlabel('$e$ [1/s]')

%% Extrapolate
% eps_t_er = linspace(-0.01, 0.06, nl);     eps_t_ext = zeros(nt, nl);                       
% for i = 1:nl, eps_t_ext(:, i) = eps_t_er(i)*ones(nt, 1);   end
% 
% l_ext = NaN(nt, nl);     gamma_s_ext = NaN(nt, nl);  
% % compute corresponding l_ext
% for j = 1:nl
%     for i = 1:nt
%         [ l_bar, ~, ~, ~, ~ ] = comp_l_bar_pol_2dof([0 theta(i)]');  
%         l_ext(i, j) = (l_bar(2) - L_t0(2)*eps_t_ext(i, j) - L_t0(2)); % needs to be adjusted for pennated muscles, 
%         % namely probably define l and then cut off depending on eps_t
% 
%         gamma_s(i, j) = comp_gs(eps_t_ext(i, j));
%     end
% end
% %
% PHI_ext     = cell(nq, 1);      LM_DOT_ext      = cell(nq, 1);   
% LN_DOT_ext  = cell(nq, 1);      LN_DOT_HAT_ext  = cell(nq, 1);
% %
% for i = 1:nq
%     for m = 1:nm
%         [Phi, ln_dot] = comp_vphi_q(nl, nt, theta, l_ext, q(i));
% 
%         Ln_dot = squeeze(ln_dot(:, :, 2))./(p8*L_0(2));     
% 
%         W_ln_dot = W{i};        ln_dot_v_hat = W_ln_dot*Phi;       
%     %     Ln_dot_hat = reshape(ln_dot_v_hat', nt, nl);
% 
%         LM_DOT_ext{i} = ln_dot;     PHI_ext{i} = Phi;       LN_DOT_ext{i} = Ln_dot;
%         LN_DOT_HAT_ext{i} = ln_dot_v_hat;
%     end
% end
% 
% eps_t_ext_v = eps_t_ext(:);
% 
% figure(31) % predicted vs actual extrapolation - deemed reasonable 16.06.2022
% for i = 1:nq
% ln_dot_v_hat = LN_DOT_HAT{i};   Ln_dot = LN_DOT_ext{i};
% subplot(4, 5, i),   surf(repmat(theta', 1, nl)/d2r, eps_t_ext, Ln_dot, 'EdgeAlpha', 0, 'FaceAlpha',0.8), grid on, hold on
%                     scatter3(theta_v(1:nth:end)/d2r, eps_t_ext_v(1:nth:end), ln_dot_v_hat(1:nth:end), 'red','filled'), hold off
% end
% xlabel('$\theta$ [deg]'), ylabel('$\epsilon_{\rm t}$ [-]'), zlabel('$\dot l$ [m/s]')
