clear all

r_e = load('data_files/r_osim.mat').r_elbow_flex;   params

n = 200;

theta     = linspace(-deg2rad(0), 2.5, n);
R         = zeros(7,n);
dr_dtheta = zeros(7,n);

p   = 12;
Phi = NaN(p+1,n);

for i=1:n
    for j = 1:p+1,    Phi(j,i) = theta(i)^(j-1);         end
end

for m = 1:nm  
    muscle_name = muscle_names(m);
    if m == 1,          R(m, :) = zeros(200, 1);
    else,               R(m,:) = r_e.(string(muscle_name));
    end
end

W7 = R*pinv(Phi); 

figure('Position', [100 100 1920 1080])
for i = 1:7
    subplot(7,1,i), plot(rad2deg(theta), R(i,:), 'b', rad2deg(theta), W7(i,:)*Phi, 'r--', 'LineWidth', 1.3), grid on, ylabel('$r_{\rm e}$ [m]'), subtitle(string(muscle_names(i)))
end
xlabel('$\theta$ [deg]'), legend('OpenSim', 'pol fit')

% save('data_files/e_r_pol','W7');
% 
% trilong_file = fopen('~/Desktop/r_models/TRIlong.dat','w');
% trilat_file = fopen('~/Desktop/r_models/TRIlat.dat', 'w');
% trimed_file = fopen('~/Desktop/r_models/TRImed.dat', 'w');
% biclong_file = fopen('~/Desktop/r_models/BIClong.dat','w');
% bicshort_file = fopen('~/Desktop/r_models/BICshort.dat', 'w');
% bra_file = fopen('~/Desktop/r_models/BRA.dat', 'w');
% 
% fprintf(trilong_file, '%f \n', W6(1, :));
% fprintf(trilat_file, '%f \n', W6(2, :));
% fprintf(trimed_file, '%f \n', W6(3, :));
% fprintf(biclong_file, '%f \n', W6(4, :));
% fprintf(bicshort_file, '%f \n', W6(5, :));
% fprintf(bra_file, '%f \n', W6(6, :));
% 
% fclose(trilong_file);
% fclose(trilat_file);
% fclose(trimed_file);
% fclose(biclong_file);
% fclose(bicshort_file);
% fclose(bra_file);


