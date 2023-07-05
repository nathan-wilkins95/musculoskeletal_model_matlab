clear all, close all

params

depst1 = 1e-3; depst2 = 1e-4; depst3 = 1e-2;
m_epst1 = -1.0:depst1:0;  m_epst2 = 1e-5:depst2:0.0201; m_epst3 = 0.0201+1e-5:depst3:0.6;
m_epst = [m_epst1 m_epst2 m_epst3];
[a, b] = size(m_epst);
gamma_s = NaN(a, b); dgs = NaN(a, b);       d2gs = NaN(a, b);

for i = 1:b
    [ gamma_s(i), dgs(i), d2gs(i) ] = comp_gs(m_epst(i));
end
    
p   = 12;
Phi = NaN(p+1,b); d_phi_d_epst = zeros(p+1,b); d2_phi_d2_epst = zeros(p+1,b);

for i=1:b
    for j = 1:p+1,    Phi(j,i) = m_epst(i)^(j-1);         end
    for j = 2:p+1,    d_phi_d_epst(j, i)   = (j-1)*m_epst(i)^(j-2);         end
    for j = 3:p+1,    d2_phi_d2_epst(j, i) = (j-1)*(j-2)*m_epst(i)^(j-3);   end
end

W = gamma_s*((Phi' * Phi + 1e-11*eye(b))\Phi');
% W = gamma_s*pinv(Phi);


figure('Position', [100 100 1920 1080])
% for i=1:7
%     subplot(7,1,i), plot(m_epst, dgs, 'b', m_epst, W*Phi, 'r--', 'LineWidth', 1.3), grid on, ylabel('$\gamma_{\rm s}$ [-]')
% end
% plot(m_epst, dgs, 'b', m_epst, W*Phi, 'r--', 'LineWidth', 1.3), grid on, ylabel('$\gamma_{\rm s}$ [-]')
subplot(311), plot(m_epst, gamma_s, 'b', m_epst, W*Phi, 'r:', 'LineWidth', 1.3), grid on, ylabel('$\gamma_{\rm s}$ [-]')
subplot(312), plot(m_epst, dgs, 'b', m_epst, W*d_phi_d_epst, 'r--', 'LineWidth', 1.3), grid on, ylabel('$\partial \gamma_{\rm s}/ \partial \epsilon_{\rm t}$ [-]')
subplot(313), plot(m_epst, d2gs, 'b', m_epst, W*d2_phi_d2_epst, 'r--', 'LineWidth', 1.3), grid on, ylabel('$\partial^2 \gamma_{\rm s}/ \partial^2 \epsilon_{\rm t}$ [-]')
xlabel('$\epsilon_{\rm t}$ [deg]')
% 

save('data_files/gs_pol','W');
