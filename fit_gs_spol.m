clear all, close all

params

eps_t_min = -0.1;       eps_t_max = 0.15;

depst1 = 1e-4; depst2 = 1e-5; depst3 = 1e-4;
m_epst1 = eps_t_min:depst1:0;  m_epst2 = 1e-5:depst2:0.0201; m_epst3 = 0.0201+1e-5:depst3:eps_t_max;
m_epst = [m_epst1 m_epst2 m_epst3];
[a, b] = size(m_epst);
gamma_s = zeros(a, b);    dgs = zeros(a, b);            d2gs = zeros(a, b);

% compute sigmoid approximation
y = zeros(a, b);          y_int = zeros(a, b);        y_dot = zeros(a, b);
y_intt = zeros(a, b);
% y_dot_int = zeros(a,b);
for i = 2:b
    [ gamma_s(i), dgs(i), d2gs(i) ] = comp_gs(m_epst(i));
    y(i) = tanh_fun(m_epst(i));
    y_intt(i) = tanh_int(m_epst(i));
    y_int(i) = y_int(i-1) + integral(@tanh_fun, m_epst(i-1), m_epst(i));
    y_dot(i) = tanh_dfun(m_epst(i));%y(i)*(51.8749 - y(i)); % 1 - y(i)^2;
%     y_dot_int(i) = y_dot_int(i-1) + 8e-5*y_dot(i);
end

% figure, plot(m_epst, y), grid on
% plot sigmoid approximation
figure('Position', [100 100 1920 1080])

subplot(411), plot(m_epst, gamma_s, 'b', m_epst, y_int, 'r--', m_epst, y_intt, 'g--', 'LineWidth', 1.3), grid on, ylabel('$\gamma_{\rm s}$ [-]'), legend('Thelen2003', 'tanh', 'Location', 'northwest')
subplot(412), plot(m_epst, dgs, 'b', m_epst, y, 'r--', 'LineWidth', 1.3), grid on, ylabel('$\partial \gamma_{\rm s}/ \partial \epsilon_{\rm t}$ [-]')
subplot(413), plot(m_epst, d2gs, 'b', m_epst, y_dot, 'r--', 'LineWidth', 1.3), grid on, ylabel('$\partial^2 \gamma_{\rm s}/ \partial^2 \epsilon_{\rm t}$ [-]')
subplot(414), plot(m_epst, y_dot, 'r--', 'LineWidth', 1.3), grid on, ylabel('$\partial^2 \gamma_{\rm s}/ \partial^2 \epsilon_{\rm t}$ [-]')
xlabel('$\epsilon_{\rm t}$ [-]')

% exportgraphics(gcf, 'results/gs_tanh.jpg', 'Resolution', 750)

% figure, plot(m_epst, y, 'b', m_epst, y_dot_int, 'r--'), grid on
% generate data for polynomial approximation
% p   = 2;
% Phi = zeros(p+1,b); Phi_int = zeros(p+1,b); d2_phi_d2_epst = zeros(p+1,b);
% 
% phi_fun = @(x) x.^(1:p+1);     
% for i=2:b
% %     for j = 1:p+1,    Phi(j,i) = m_epst(i)^(j-1);         end
% % %     for j = 2:p+1,    d_phi_d_epst(j, i)   = (j-1)*m_epst(i)^(j-2);         end
% % %     for j = 3:p+1,    d2_phi_d2_epst(j, i) = (j-1)*(j-2)*m_epst(i)^(j-3);   end
%     Phi(:, i)               = phi_fun(m_epst(i));
%     Phi_int(:, i)         = Phi_int(:, i-1) + integral(phi_fun, m_epst(i-1), m_epst(i), 'ArrayValued', true)';
% end
% 
% PHI = [Phi_int, Phi];
% W = ([gamma_s, dgs]' - [y_int, y]')' * pinv(PHI); % ((PHI' * PHI + 1e-10*eye(2*b))\PHI');
% 
% GS_SIG_hat = W*PHI;
% 
% 
% figure('Position', [100 100 1920 1080])
% 
% subplot(411), plot(m_epst, gamma_s, 'b', m_epst, y_int, 'r--', m_epst, GS_SIG_hat(1:b) + y_int, 'g-.', 'LineWidth', 1.4), grid on, ylabel('$\gamma_{\rm s}$ [-]'), legend('Thelen2003', 'SIG', 'SIG+POl', 'Location', 'northwest')
% subplot(412), plot(m_epst, dgs, 'b', m_epst, y, 'r--', m_epst, GS_SIG_hat(b+1:end) + y, 'g-.', 'LineWidth', 1.4), grid on, ylabel('$\partial \gamma_{\rm s}/ \partial \epsilon_{\rm t}$ [-]')
% subplot(413), plot(m_epst, d2gs, 'b', m_epst, y_dot, 'r--', 'LineWidth', 1.3), grid on, ylabel('$\partial^2 \gamma_{\rm s}/ \partial^2 \epsilon_{\rm t}$ [-]')
% subplot(414), plot(m_epst, y_dot, 'r--', 'LineWidth', 1.3), grid on, ylabel('$\partial^2 \gamma_{\rm s}/ \partial^2 \epsilon_{\rm t}$ [-]')
% xlabel('$\epsilon_{\rm t}$ [-]')

% exportgraphics(gcf, 'results/gs_sig_pol.jpg', 'Resolution', 1000)


% syms f(x)
% a = 300;        c = 51.8749;         d = 0.0135;     a1 = 0.25;   
% f(x) = c*( exp(a*(x - d)) )./( exp(a*(x - d)) + exp(-a1*a*(x - d)) );
% 
% F = int(f, x)
% df = diff(f, x)
