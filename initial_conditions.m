disp('Initialisation Start')

% Initial Joint Angle, Velocity and Torque
% theta   = theta_des+0.01;      theta_dot = theta_dot_des;      
tau0    = [0 0]';               
theta   = [0 0]';         theta_dot = [0 0]';

[ l_bar0, dl_bar_dtheta0, ~, ~, ~ ] = comp_l_bar_pol_2dof(theta);          
[R0, ~, ~]   = comp_r_pol(theta); 
Theta0       = [theta' theta_dot']';
%% 
l       = L_0;
eps_t   = zeros(nm,1);  gs = zeros(nm,1);  lm_dot = zeros(nm, 1);

m_l_t   = L_t0*( 1 + 0.016); %0.0105                                       % stretched tendons initially
m_l     = l_bar0  - m_l_t;

for i = 1:nm
    eps_t(i)   = ( m_l_t(i) - L_t0(i)  ) / L_t0(i);
    gs(i) = comp_gs_tpol(eps_t(i));
end
%
tau_des0  = [5 10]';
%
%% Find eps_t to gamma_sdes
n =  101;        tw_init = linspace(0,1e-1,n);

[ t_init, eps_t_init ] = ode45(@(t,y) trim_dyn(t, y, theta, tau_des0), tw_init, eps_t, options);   

tau_init = zeros(2, n);
eps_t0   = eps_t_init(end,:)';

for i = 1:n
    gs = zeros(nm,1);
    for j = 1:nm, gs(j) = comp_gs_tpol(eps_t_init(i,j));         end
    tau_init(:, i) = R0'*Fmax*gs;
end

% figure, plot(t_init,tau_init), grid on, xlabel('$t$ [s]'), ylabel('$\tau(t)$ [N.m]'), legend('\tau_{\rm e}', '\tau_{\rm s}')
% % exportgraphics(gcf, savedir + "tau.png", 'Resolution', res)
% 
% figure
% for i=1:nm
%     subplot(nm,1,i), plot(t_init,eps_t_init(:,i),'b',t_init([1,end]),eps_th*[1,1],'g--',t_init([1,end]),0*[1,1],'g--'), grid on, ylabel('$\epsilon_{\rm  t}$ [-]'),
% end
% subplot(nm,1,nm), xlabel('$t$ [s]')
% % exportgraphics(gcf, savedir + "eps_t.png", 'Resolution', res)

%% Find muscle length to complete mt length
lt = diag(L_t0)*( ones(nm,1) + eps_t0 );
l0 = L_0;

[ t_i2, l_i ] = ode45(@(t,y) trim_length(t, y, l_bar0, lt), tw_init, l0, options);

lca = zeros(n, nm);
for i = 1:n
    c_alpha  = zeros(nm,1);
    for j = 1:nm,    c_alpha(j) = sqrt( 1 - ( L_0(j)*sin(Alpha_0(j))/l_i(i,j) )^2 );    end
    lca(i,:) = ( diag(l_i(i,:))*c_alpha )';
end

% figure
% for i=1:nm
%     subplot(nm,1,i), plot(t_i2,lca(:,i),'b',t_i2([1,end]),(l_bar0(i) - lt(i))*[1,1],'g--'), grid on, ylabel('$l_0 \cos(\alpha)$ [m]'),
% end
% subplot(nm,1,nm), xlabel('$t$ [s]')
% 
% figure
% for i=1:nm
%     subplot(nm,1,i), plot(t_i2,l_i(:,i),'b',t_i2([1,end]),L_0(i)*[1,1],'g--'), grid on, ylabel('$l_0$ [m]'),
% end
% subplot(nm,1,nm), xlabel('$t$ [s]')
% 
% figure
% for i=1:nm
%     subplot(nm,1,i), plot(t_i2,lca(:,i)+lt(i),'b',t_i2([1,end]),l_bar0(i)*[1,1],'g--'), grid on, ylabel('$l_0+l_{{\rm t}0}$ [m]'),
% end
% subplot(nm,1,nm), xlabel('$t$ [s]')

%% Find activation to balance muscle length dynamics
l0          = l_i(end,:);           eps_t0      = eps_t_init(end,:);
n           = 101;                 tw_init     = linspace(0,1e-3,n);
%
c_alpha0    = zeros(nm,1);
for i=1:nm,     c_alpha0(i)  = sqrt( 1 - ( L_0(i)*sin(Alpha_0(i))/l0(i) )^2 );     end
%
q0      = 0.1*ones(nm,1);
% 
[ t_i2, q_i ]  = ode45(@(t, q) trim_musc(t, q, l0, c_alpha0, eps_t0 ), tw_init, q0, options);                  % lm_dot -> 0 
% [ t_i2, q_i ]  = ode45(@(t, q) trim_musc_e4(t, q, l0, c_alpha0, eps_t0, lm_dot_des), tw_init, q0, options);    % e4 -> 0   
%%
ga = zeros(nm,1);          gp = zeros(nm,1);
gs = zeros(nm,1);          gc = zeros(nm,1);
ln_dot  = zeros(n, nm);

for j = 1:nm
        l_n                        = l0(j) / L_0(j);
        [ gp(j), ga(j) ] = comp_g_pa(l_n);
        gs(j)                 = comp_gs_tpol(eps_t0(j));
        gc(j)                 = gs(j)/c_alpha0(j) - gp(j);
    for i = 1:n
        ln_dot(i,j) = ncomp_ln_dot([ga(j), gc(j)], q_i(i,j));
    end
end

% figure
% for i=1:nm
%     subplot(nm,1,i), plot(t_i2,q_i(:,i),'b',t_i2([1,end]),0.01*[1,1],'g--',t_i2([1,end]),[1,1],'g--'), grid on, ylabel('$q_{0}$ [-]'),
% end
% subplot(nm,1,nm), xlabel('$t$ [s]')
% 
% figure
% for i=1:nm
%     subplot(nm,1,i), plot(t_i2,ln_dot(:,i),'b',t_i2([1,end]),0*[1,1],'g--'), grid on, ylabel('$\dot l_{{\rm n}0}$ [1/s]'),
% end
% subplot(nm,1,nm), xlabel('$t$ [s]')
q0 = q_i(end,:)';

disp('Initialisation Done')

% Recomp for coherence
[ ln_dot0, Gamma, meps_t, c_alpha, ~, mpartials, ~ ] = comp_mus_2dof( theta, l0, q0 );
%
ln_dot0 = ln_dot0';     
gs      = Gamma(:,3);
%