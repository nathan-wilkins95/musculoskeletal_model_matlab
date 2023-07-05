%%
global eps_th a0 a1 a2 a3 p1 p2 p3 p4 p5 p6 p7 p8 k e_0m g ...
       tau_a tau_d alpha_1 m_specs muscle_names l1 l2 m1 m2 gr ...
       W7_r W6_l_bar Wp Ws Wc p nm nj Fmax L_0 L_t0 Alpha_0 ...
       Sigma Sigma_l Sigma_q Sigma_f des_params d2r ...
       eps_p q_exp rho sigma q_min q_max dt tspan options tau_f WSc ...
       theta_sf theta_ef theta_si theta_ei T         


addpath(genpath('./ms_model/')),    addpath(genpath('./plot_utils/'))
addpath(genpath('./data_utils/')),  addpath(genpath('./utils/'))
addpath(genpath('./tests/'))

% muscle dynamics parameters
p1 = 1/4;                       p2 = 3/4;                       p3 = 1.8;
p4 = 0.95;                      p5 = 3/10;                       
p6 = p3 - 1;                    p7 = 2 + 2 / p5;                p8 = 1e1;                    

% tendon gain parameters
eps_0 = 0.033;

eps_th = 99*eps_0*exp(3) / (166*exp(3)-67);
a2     = 67 / ( 100*(eps_0 - (99*eps_0*exp(3))/(166*exp(3)-67)));
a1     = 0.33;      a3 = 3;

% muscle gain parameters
a0      = 0.5;              % 0.6 for young
k       = 4;
e_0m    = 0.6;             g  = 0.5;

% q_dot parameters
tau_a = 0.01;       tau_d = 0.04;   q_min = 0.01;   q_max = 1;      

alpha_1 = acos(0.1);                                                       % maximum pennation angle

tau_f = 0.001;
% muscle specific parameters
muscle_data = load('data_files/m_specs.mat');
m_specs = muscle_data;      muscle_names = fieldnames(m_specs);

nm   = length(muscle_names);      nj      = 2;
Fmax = zeros(nm,1);               L_0     = zeros(nm,1);
L_t0 = zeros(nm,1);               Alpha_0 = zeros(nm,1);

for i = 1:nm
    muscle_name = muscle_names{i};
    Fmax(i)     = m_specs.(string(muscle_name)).params.fmax;               % maximum isometric force 
    L_0(i)      = m_specs.(string(muscle_name)).params.opt_fiber_length;   % optimal fiber length
    L_t0(i)     = m_specs.(string(muscle_name)).params.lt_slack;           % tendon slack length
    Alpha_0(i)  = m_specs.(string(muscle_name)).params.alpha;              % pennation angle at optimal 
end
%
Fmax = diag(Fmax); % maximum isometric force
% polynomial coefficients
w7_e = load('data_files/e_r_pol.mat');              W7_e = w7_e.W7;        % elbow joint
w7_s = load('data_files/s_r_pol.mat');              W7_s = w7_s.W7;        % shoulder joint
w6_l_bar = load('data_files/l_bar_pol_2dof.mat');   W6_l_bar = w6_l_bar.W6;% l_bar
wp = load('data_files/gamma_pol.mat');              Wp = wp.Wp;            % gamma_p
wc = load('data_files/c_alpha_pol.mat');            Wc = wc.W;             % c_alpha
ws = load('data_files/gs_pol');                     Ws = ws.W;             % gamma_s
p = 5;                                              W7_r = {W7_s W7_e};


% Sigma 0.25
Sigma   = 1e-2*eye(nm);     Sigma_l = 1e4*eye(nm);       Sigma_q = 1e5*eye(nm);
Sigma_f = 1e5*eye(nm);                                                     % step sizes for gradient descent (initial conditions)
   
% skeletal parameters 
m1 = 1.86;      m2 = 1.53;  % segment weight [kg]
l1 = 0.29;      l2 = 0.24;  % segment length [m]

gr = 9.81;      % gravity [m/s^2]
d2r = pi / 180; % rad/deg conversion parameter

% projection parameters
eps_p = 0.5;   q_exp = 6;  rho = 0.5;  sigma = 0.49;

% desired trajectory parameters
a_e_des     = 45*d2r;               a_s_des = 0.25*a_e_des;
b_e_des     = 1.5*a_e_des;          b_s_des     = 0.5*a_e_des;
w_e_des     = 8/(2*pi);             w_s_des = 6/(2*pi);
phi_e_des   = -4*(pi/2);            phi_s_des   = -4*(pi/2);
%             
des_params = [ a_e_des b_e_des w_e_des phi_e_des a_s_des b_s_des w_s_des phi_s_des];

WSc = [1        -0.5    0       0       0       0       0
       -0.5     1       0       0       -0.5    0       0
       0        0       1       0       0       -0.5    -0.5
       0        0       0       1       0       -0.5    -0.5
       0        -0.5    0       0       1       0       0
       0        0       -0.5    -0.5    0       1       0
       0        0       -0.5    -0.5    0       0       1];


% flexion trajectory parameters
theta_sf = 45;  theta_si = 0;  theta_ef = 85; theta_ei = 0; T = 2; 

% matlab settings
set(0,'defaulttextInterpreter','latex')
options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06); 

dt = 0.001;     tspan =(0:dt:10);