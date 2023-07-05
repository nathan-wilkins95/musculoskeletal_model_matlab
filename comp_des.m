function [theta_des, theta_des_dot, theta_des_ddot, theta_des_3dot, theta_des_4dot] = comp_des(t)

global des_params

%shoulder params
a_s = des_params(5);        b_s = des_params(6);
w_s = des_params(7);        phi_s = des_params(8);

% elbow params
a_e = des_params(1);        b_e   = des_params(2);
w_e = des_params(3);        phi_e = des_params(4);

a = [a_s a_e]';     w = [w_s w_e]';     b = [b_s b_e]';     phi = [phi_s phi_e]';


theta_des       = a.*sin(t*w + phi) + b;
theta_des_dot   = a.*w.*cos( t*w + phi );
theta_des_ddot  = -a.*(w.^2).*sin( t*w + phi );
theta_des_3dot  = -a.*(w.^3).*cos( t*w + phi );
theta_des_4dot  = a.*(w.^4).*sin( t*w + phi );

% theta1_des      =         a_s*sin( w_s*t + phi_s ) + b_s;
% theta1_des_dot  =       a_s*w_s*cos( w_s*t + phi_s );
% theta1_des_ddot =  -a_s*(w_s^2)*sin( w_s*t + phi_s );
% theta1_des_3dot =  -a_s*(w_s^3)*cos( w_s*t + phi_s );
% theta1_des_4dot =   a_s*(w_s^4)*sin( w_s*t + phi_s );
% 
% theta2_des      =         a_e*sin( w_e*t + phi_e )  + b_e;
% theta2_des_dot  =       a_e*w_e*cos( w_e*t + phi_e );
% theta2_des_ddot =  -a_e*(w_e^2)*sin( w_e*t + phi_e );
% theta2_des_3dot =  -a_e*(w_e^3)*cos( w_e*t + phi_e );
% theta2_des_4dot =   a_e*(w_e^4)*sin( w_e*t + phi_e );
% 
% % theta_des           = [theta1_des theta2_des]';
% theta_des_dot       = [theta1_des_dot theta2_des_dot]';
% theta_des_ddot      = [theta1_des_ddot theta2_des_ddot]';
% theta_des_3dot      = [theta1_des_3dot theta2_des_3dot]';
% theta_des_4dot      = [theta1_des_4dot theta2_des_4dot]';

end