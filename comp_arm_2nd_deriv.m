function g_ddot = comp_arm_2nd_deriv(theta, theta_dot, theta_ddot, theta_3dot)
% Compute second order time derivatives of skeletal dynamics - verified

global m1 m2 l1 l2 gr

t1 = theta(1);                  t2 = theta(2);
w1 = theta_dot(1);              w2 = theta_dot(2);
theta_s_ddot = theta_ddot(1);   theta_e_ddot = theta_ddot(2);
theta_s_3dot = theta_3dot(1);   theta_e_3dot = theta_3dot(2);
%
C       = [0 (l2 / l1)*(m2/(m1+m2))*sin(t1-t2)
            -(l1/l2)*sin(t1-t2) 0];  

dC1_dt1 = (l2 / l1) * (m2 / (m2 + m1))*cos(t1-t2)*w1;       dC1_dt2 = -(l2 / l1) * (m2 / (m2 + m1))*cos(t1-t2)*w2;
dC2_dt1 = -(l1/l2)*cos(t1-t2)*w1;                           dC2_dt2 = (l1/l2)*cos(t1-t2)*w2;

C_dot   = [0 dC1_dt1+dC1_dt2
           dC2_dt1+dC2_dt2  0];
%
dC1_2_dt1_2     = -(l2 / l1) * (m2 / (m2 + m1))*(sin(t1-t2)*w1^2 - cos(t1-t2)*theta_s_ddot);  
dC1_2_dt2_2     = -(l2 / l1) * (m2 / (m2 + m1))*(sin(t1-t2)*w2^2 + cos(t1-t2)*theta_e_ddot);
dC1_2_dt1dt2    = (l2 / l1) * (m2 / (m2 + m1))*sin(t1-t2)*w1*w2;

dC2_2_dt1_2     = (l1/l2)*(sin(t1-t2)*w1^2 - cos(t1-t2)*theta_s_ddot);       
dC2_2_dt2_2     = (l1/l2)*(sin(t1-t2)*w2^2 + cos(t1-t2)*theta_e_ddot);
dC2_2_dt1dt2    = -(l1/l2)*sin(t1-t2)*w1*w2;

C_ddot  = [0 (dC1_2_dt1_2 + dC1_2_dt2_2 + 2*dC1_2_dt1dt2)
            (dC2_2_dt1_2 + dC2_2_dt2_2 + 2*dC2_2_dt1dt2) 0];
%       
%     gt      = [gr/l1*sin(t1) gr/l2*sin(t2)]';                               
%     gt_dot  = [(gr/l1)*cos(t1)*w1 (gr/l2)*cos(t2)*w2]';
gt_ddot = [-(gr/l1)*(sin(t1)*w1^2 - cos(t1)*theta_s_ddot), -(gr/l2)*(sin(t2)*w2^2 - cos(t2)*theta_e_ddot)]';
%
%     g_dot = C_dot * [w1^2 w2^2]' + C * [2*w1*theta_s_ddot 2*w2*theta_e_ddot]' + gt_dot ;
g_ddot =    C_ddot*[w1^2 w2^2]' + 4*C_dot*[w1*theta_s_ddot w2*theta_e_ddot]' + ...
            2*C*([theta_s_ddot^2 theta_e_ddot^2]' + [w1*theta_s_3dot w2*theta_e_3dot]' ) + gt_ddot;


end