clear all; close all;clc
% file to test projection algorithm (Pomet, Praly)
global eps_p q_exp rho sigma

% projection parameters
eps_p = 0.5;   q_exp = 8;  rho = 0.5;  sigma = 0.49;      

dt = 0.001;     nm = 7;
tspan =(0:dt:10);

% matlab settings
set(0,'defaulttextInterpreter','latex')
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10); 


q0 = 0.5;   q_dot0 = 1e2;   % throw q at boundary and see if proj prevents it from crossing the boundary
X0 = [q0 q_dot0]';

[t, X] = ode45(@dynamics_test_phi, tspan, X0, options);

nn = length(X); % phi_dot = NaN(nn, 1);   Phi = NaN(nn, 1);

figure(1)   % check to see that despite q_dot increasing, q does not cross boundary
plot(t, X(:, 1)), xlabel('$time$ [s]'), grid on, ylabel('$q$ [-]')   

figure(2)
plot(t, X(:, 2)), xlabel('$time$ [s]'), grid on, ylabel('$\dot q$ [-]')
