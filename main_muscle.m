clear all, close all, clc

global l_0 l_t0 alpha_0
params
% BRA params to make main.m work
l_0     = 0.0858;
l_t0    = 0.0535;
alpha_0 = 0.0;
%% Define range of inputs
mpc      = 2/100;
l_r      = [ 1-mpc 1+mpc ]*l_0;                             nl     = 101;
l_bar_r  = [ 1-mpc 1+mpc ]*(l_0 + l_t0);                  nl_bar = 101;
q_r      = [ 0.01 1 ];                                         nq     = 101;

l_vec     = minterp(l_r,nl);
l_bar_vec = minterp(l_bar_r,nl_bar);
q_vec     = minterp(q_r,nq);

%% Generate sample surface

Gamma_p = NaN(nl,nl_bar);         Gamma_a = NaN(nl,nl_bar);
Gamma_s = NaN(nl,nl_bar);         Gamma_c = NaN(nl,nl_bar);
ln_dot  = NaN(nl,nl_bar);         Mflag   = NaN(nl,nl_bar);
Eps_t   = NaN(nl,nl_bar);         Calpha  = NaN(nl,nl_bar);
q       = 1;

for i = 1:nl
    for j = 1:nl_bar
            x = [ l_vec(i) l_bar_vec(j) q ];
            
            % normalized muscle length
            l_n = l_vec(i) / l_0;

            % passive & active force gains
            [ gamma_p, gamma_a ] = comp_g_pa(l_n);

            % series force gain
            c_alpha = sqrt( 1 - ( l_0*sin(alpha_0)/l_vec(i) )^2 );
            Eps_t(j, i)   = ( l_bar_vec(j) - l_vec(i)*c_alpha - l_t0) / l_t0;
            gamma_s = comp_gs(Eps_t(j, i));

            % contractile force gain
            gamma_c = gamma_s/c_alpha - gamma_p;

            ln_dot(j,i) = comp_ln_dot(gamma_c, gamma_a, q);
            Gamma_p(j,i) = gamma_p;         Gamma_a(j,i) = gamma_a;
            Gamma_s(j,i) = gamma_s;         Gamma_c(j,i) = gamma_s;
            Calpha(j,i)  = c_alpha;
    end
end

%% Check gamma_s
ne = nl*1000;
eps_t_vec   = linspace(0,eps_0*1.1,ne);
gamma_s_vec = NaN(nl,1);

for i=1:ne
    gamma_s_vec(i) = comp_gs(eps_t_vec(i));
end

figure('Position', [100 100 400 300]), plot(eps_t_vec,gamma_s_vec, 'LineWidth', 1.4), hold on, plot([ eps_t_vec(1), eps_th ],a1*ones(1,2),'r:'), hold on
           plot(eps_th*ones(1,2),[0,a1],'r:',eps_th*ones(1,2),[0,a1],'r.'), hold on
           plot(eps_0*ones(1,2),[0,1],'r:',[0,eps_0],[1,1],'r:',eps_0*ones(1,2),[0,1],'r.')
           grid on, xlabel('$\epsilon_{\rm t}$ [-]', 'FontSize', 22), ylabel('$\gamma_{\rm s}$ [-]', 'FontSize', 22)%, title('Tendon gain as a function of strain'), hold off 
%saveas(gcf, './results/gamma_s.png')
%% Check ln_dot
ngc        = ne;
gc_vec     = linspace(-0.1,p3*1.1,ngc);
ln_dot_vec = NaN(nl,1);             mf     = NaN(nl,1);
mga        = 1;                     mq     = 1;

for i=1:ngc
    [ ln_dot_vec(i), mf(i) ] = comp_ln_dot(gc_vec(i),mga,mq);
end

figure(2), plot(gc_vec,ln_dot_vec, 'LineWidth', 1.3), hold on, plot([ gc_vec(1), p3 ],ones(1,2),'r--',p3,1,'r.'), hold on
           plot([ gc_vec(1), p3 ],zeros(1,2),'r--',p3,0,'r.'), hold on
           plot(ones(1,2),[-1,0],'r:',1,0,'r.'), hold on
           plot(p3*ones(1,2),[-1,1],'r--',p3,1,'r.'), hold on
           plot(mq*mga*p3*p4*ones(2,1),[-1,1],'k--',mq*mga*ones(2,1),[-1,1],'k--')       
           grid on, xlabel('$\gamma_{\rm c}$ [-]', 'FontSize', 22), ylabel('$\dot l_{\rm n}$ [1/s]', 'FontSize', 22), ylim([-1.0, 1.25]), xlim([0,1.85])%, title('Velocity as a function of muscle force gain')
           hold off 
%saveas(gcf, './results/ln_dot.png')
%% Check Gamma_a, gamma_p
mn       = 1001;
m_ln     = linspace(0,2,mn);
mGamma_p = NaN(1,mn);           mGamma_a = NaN(1,mn);

for i = 1:mn
    [ mGamma_p(i), mGamma_a(i) ] = comp_g_pa(m_ln(i));
end

figure('Position', [200 200 500 400]), plot(m_ln,mGamma_p,'k',m_ln,mGamma_a,'b', 'LineWidth', 1.4), hold on
           plot([m_ln(1),(1+e_0m)],ones(1,2),':k',(1+e_0m)*ones(2,1),[0,1],'k:',1+e_0m,1,'r.')
           grid on, xlabel('$l_{\rm n}$ [-]', 'FontSize', 22), ylabel('$\gamma_{\rm p}$, $\gamma_{\rm a}$ [-]', 'FontSize', 22), xlim([0.2,1.8]), ylim([0 1.5])%, title('Muscle force gains as a function of length')
           hold off
           %saveas(gcf, './results/gamma_a_p.png')
%% Display
% figure, surf(l_vec,l_bar_vec,Eps_t), grid on, xlabel('muscle length $l$ [m]'), ylabel('total length $\bar l$ [m]'), zlabel('tendon strain $\epsilon_{\rm t}$ [-]'), %zlim ([-0.01,eps_0+0.01]) 

% figure, contourf(l_vec,l_bar_vec,Eps_t,[0,eps_th,eps_0],'showtext','on'), hold on,
%         plot([ l_vec(1), l_vec(end)], (l_0+l_t0)*ones(1,2),'r--'), hold on
%         plot(l_0*ones(1,2), [ l_bar_vec(1), l_bar_vec(end)],'r--'), xlabel('muscle length $l$ [m]'), ylabel('total length $\bar l$ [m]'), hold off

figure, surf(l_vec,l_bar_vec,Gamma_s), grid on, xlabel('muscle length $l$ [m]'), ylabel('total length $\bar l$ [m]'), zlabel('series gain $\gamma_{\rm s}$ [-]') , %zlim ([0,1.1])

% figure, contourf(l_vec,l_bar_vec,Gamma_s,[0,a1,1],'showtext','on'), hold on,
%         plot([ l_vec(1), l_vec(end)], (l_0+l_t0)*ones(1,2),'r--'), hold on
%         plot(l_0*ones(1,2), [ l_bar_vec(1), l_bar_vec(end)],'r--'), xlabel('muscle length $l$ [m]'), ylabel('total length $\bar l$ [m]'), hold off

% figure, surf(l_vec,l_bar_vec,Gamma_c), grid on, xlabel('muscle length $l$ [m]'), ylabel('total length $\bar l$ [m]'), zlabel('contractile gain $\gamma_{\rm c}$ [-]') 

% figure, surf(l_vec,l_bar_vec,Gamma_c-Gamma_s), grid on, xlabel('$l$ [m]'), ylabel('$\bar l$ [m]'), zlabel('$\Delta \gamma$ [-]') 

figure, surf(l_vec,l_bar_vec,ln_dot), grid on, xlabel('muscle length $l$ [m]'), ylabel('total length $\bar l$ [m]'), zlabel('length RoC $\dot l_{\rm n}(t)$ [1/s]'), zlim([-1.1 1.1])
% figure, surf(l_vec,l_bar_vec,Mflag), grid on, xlabel('muscle length $l$ [m]'), ylabel('total length $\bar l$ [m]'), zlabel('flag')