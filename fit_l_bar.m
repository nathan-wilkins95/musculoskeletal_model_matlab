clear all

l_bar_struct = load('data_files/l_bar_data.mat');
muscle_data = load('data_files/m_specs.mat');
m_specs = muscle_data;      muscle_names = fieldnames(m_specs);

n = 100;

theta_s     = linspace(-pi/2, pi/2, n);     theta_e     = linspace(-deg2rad(5), 2.7, n); % as in pycharm
L_bar         = zeros(6, n*n);
[Te, Ts] = meshgrid(theta_s, theta_e);

mdls = cell(7, 1);

for m = 1:7  
    muscle_name = muscle_names(m);
    data = l_bar_struct.(string(muscle_name));
    L_bar(m, :) = data(:);
    mdls{m} = fitrgp([Ts(:) Te(:)], L_bar(m, :), 'KernelFunction', 'ardexponential');
    m
end

Ts = Ts(:); Te = Te(:);
%% Visualization

figure('Position', [100 100 1920 1080])
for i = 1:7
        mdl = mdls{i};
%         sf = mdl.KernelInformation.KernelParameters(3);
%         sl1 = mdl.KernelInformation.KernelParameters(1);
%         sl2 = mdl.KernelInformation.KernelParameters(2);
%         alpha = mdl.Alpha;  beta = mdl.Beta;
%         AS = mdl.ActiveSetVectors;
%         K = NaN(length(Ts), length(AS));
%         for k = 1:length(Ts)
%             
% %                 sqdist = sqrt(sum(([Ts(k) Te(k)] - AS(j, :)).^2) / sl^2);
%                 sqdist = sqrt(sum(([Ts(k) Te(k)] - AS).^2 ./ [sl1 sl2].^2, 2));
%                 K(k, :) = sf^2 * exp(-sqdist);
%             
%         end
%         y_hat = beta + K * alpha;
    y_hat = predict(mdl, [Ts Te]);
    subplot(2, 4, i), surf(rad2deg(theta_s), rad2deg(theta_e), reshape(L_bar(i, :)', n,n)', 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'FaceColor', 'blue'), hold on, ylabel('$\theta_{\rm e}$ [deg]'), xlabel('$\theta_{\rm s}$ [deg]'), zlabel('$\bar l$ [m]')
                      surf(rad2deg(theta_s),  rad2deg(theta_e), reshape(y_hat, n,n)', 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'FaceColor', 'red'), hold off
                                subtitle(string(muscle_names(i)))
    mean((reshape(y_hat, n, n) - reshape(L_bar(i, :), n, n)).^2, 'all')
end

exportgraphics(gcf, "~/Nextcloud/Documents/PhD/Yannick/arm_control/backstepping/l_bar_pol/l_bar_pol.png", 'Resolution', 700)


%% Save models

mp = cell(7, 1);
for i = 1:7
    
    mdl = mdls{i};      
    sf = mdl.KernelInformation.KernelParameters(3);                        % sigma_f, kernel parameter
    sl1 = mdl.KernelInformation.KernelParameters(1);                       % sigma_l1, kernel parameter, 1st predictor
    sl2 = mdl.KernelInformation.KernelParameters(2);                       % simga_l2, kernel parameter, 2nd predictor
    alpha = mdl.Alpha;  beta = mdl.Beta;                                   % weights; offset
    AS = mdl.ActiveSetVectors;                                             % active set for covariance matrix
    
    p = struct('sigma_f', sf, 'sigma_l1', sl1, 'sigma_l2', sl2, ...
        'alpha', alpha, 'AS', AS, 'beta', beta);
    mp{i} = p;
end
save('data_files/l_bar_rgp','mp');


