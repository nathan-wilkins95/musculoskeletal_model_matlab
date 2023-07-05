clear all, clc, close all

global options L_0 Fmax filt

tic, disp('working...');
params
%%
% Initial conditions
initial_conditions

X0 = [ l0 q0' Theta0']';

[t, X] = ode45(@dynamics, tspan, X0, options); toc
disp('done');

%% Recompute

theta       = X(:, 15:16);         theta_dot       = X(:, 17:18);
l_int       = X(:, 1:7);           q_int           = X(:, 8:14);

nn = length(t);

%% Plotting

figure(1), tiledlayout(2, 1,"TileSpacing","compact")
nexttile,   plot(t, theta/d2r), grid on, ylabel('$\theta$ [deg]'), legend('$\theta_{\rm s}$', '$\theta_{\rm e}$', 'Interpreter', 'latex')
nexttile,   plot(t, theta_dot/d2r), grid on, ylabel('$\dot \theta$ [deg/s]'), legend('$\theta_{\rm s}$', '$\theta_{\rm e}$', 'Interpreter', 'latex')
xlabel('time [s]')

figure(2), tiledlayout(nm, 1, "TileSpacing","compact")
for i = 1:nm
    nexttile, plot(t, l_int(:, i)), grid on, ylabel('$l$ [m]'), legend(muscle_names{i})
end, xlabel('time [s]')

figure(3), tiledlayout(nm, 1, "TileSpacing","compact")
for i = 1:nm
    nexttile, plot(t, q_int(:, i)), grid on, ylabel('$q$ [-]'), ylim([0 1]), legend(muscle_names{i})
end, xlabel('time [s]')