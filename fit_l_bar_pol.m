clear all

l_bar_struct = load('data_files/l_bar_data.mat');
muscle_data = load('data_files/m_specs.mat');
m_specs = muscle_data;      muscle_names = fieldnames(m_specs);

n = 100; nm = 7;


theta_s     = linspace(-pi/2, pi/2, n);
% theta_s     = (theta_s - min(theta_s))/(max(theta_s) - min(theta_s))*2 - 1;
theta_e     = linspace(-deg2rad(5), 2.7, n);
% theta_e     = (theta_e - min(theta_e))/(max(theta_e) - min(theta_s));
L_bar         = zeros(nm, n*n);

p_vec = 0:1:7;
Phi = [];

for m = 1:nm 
    muscle_name = muscle_names(m);
    data = l_bar_struct.(string(muscle_name));
    L_bar(m, :) = reshape(data.', 1, []); % make sure to flatten by rows
end

for i = 1:n 
    for j = 1:n
        Phi = [ Phi comp_vphi(theta_s(i), theta_e(j), p_vec)' ];
    end 
end


W6 = L_bar * ((Phi' * Phi + 1*eye(10000))\Phi');


[Ts, Te] = meshgrid(theta_s, theta_e);
for i = 1:7
    e = W6(i, :)*Phi - L_bar(i, :);

    figure('Position', [100 100 1920 1080])
    subplot(1,2,1),     scatter3(rad2deg(Te(:)), rad2deg(Ts(:)), W6(i, :)*Phi, 25, 'r.'), hold on, grid on
                        surf(rad2deg(Te), rad2deg(Ts), reshape(L_bar(i, :), n, n)), alpha 0.95, hold off
                        xlabel('$\theta_{\rm e}$ [deg]'), ylabel('$\theta_{\rm s}$ [deg]'), zlabel('$\bar l$ [m]')
                        title(muscle_names{i})
                        
    subplot(1,2,2),     surf(rad2deg(Te), rad2deg(Ts), reshape(e, n, n)), grid on 
                        xlabel('$\theta_{\rm e}$ [deg]'), ylabel('$\theta_{\rm s}$ [deg]'), zlabel('$\delta \bar l$ [m]')

%     mean((reshape(L_bar(i, :), n, n) - reshape(W6(i, :)*(Phi), n, n)).^2, 'all')
%     max(max(abs(e)))
%     exportgraphics(gcf, "~/Nextcloud/Documents/PhD/Yannick/arm_control/backstepping/l_bar_pol/"+muscle_names{i+1}+".png", 'Resolution', 600)
end

% save('data_files/l_bar_pol_2dof','W6'); 