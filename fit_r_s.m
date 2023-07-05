clear all

r_s = load('data_files/r_osim.mat').r_shoulder_elev;   params

n = 200;

theta     = linspace(-pi/2, pi/2, n)';
R         = zeros(6,n);
dr_dtheta = zeros(6,n);

p   = 12;
Phi = NaN(p+1,n);

for i=1:n
    for j = 1:p+1,    Phi(j,i) = theta(i)^(j-1);         end
end

for m = 1:7
    muscle_name = muscle_names(m);
    if (m == 3) || (m == 4) || (m == 7)                           % setting muscles not acting on shoulder to zero
        R(m, :) = zeros(200, 1);
    else
        R(m,:) = r_s.(string(muscle_name));
    end
    
    if m == 5
        mdl = fitrgp(theta, R(m, :), 'KernelFunction', 'exponential');
    end

end

W7 = R*pinv(Phi); 

figure('Position', [100 100 1920 1080])
for i = 1:7
        subplot(7,1,i), plot(rad2deg(theta), R(i,:), 'b', rad2deg(theta), W7(i, :)*Phi, 'r--', 'LineWidth', 1.3), grid on, ylabel('$r_{\rm s}$[m]'), subtitle(string(muscle_names(i)))
        if i == 5, hold on, plot(rad2deg(theta), predict(mdl, theta)), end
end
xlabel('$\theta$ [deg]'), legend('OpenSim', 'pol fit')
exportgraphics(gcf, "~/Nextcloud/Documents/PhD/Yannick/arm_control/backstepping/r_s/r_s.png", 'Resolution', 700)

save('data_files/s_r_pol','W7');

save('mdl.mat', 'mdl')