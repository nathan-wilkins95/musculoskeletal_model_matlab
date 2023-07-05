function M_ddot = comp_M_ddot(theta, theta_dot, theta_ddot)
% compute second time derivative of inertia matrix M
global m1 m2 l1 l2

t1 = theta(1);                  t2 = theta(2);
w1 = theta_dot(1);              w2 = theta_dot(2);
theta_s_ddot = theta_ddot(1);   theta_e_ddot = theta_ddot(2);

% first order time derivatives for consistency
%     a1 = (l2 / l1) * (m2 / (m1 + m2)) * cos(t1 - t2);
%     da1_dt1 = - (l2 / l1) * (m2 / (m1 + m2)) * sin(t1 - t2) * w1; da1_dt2 = (l2 / l1) * (m2 / (m1 + m2)) * sin(t1 - t2) * w2;
%     a1_dot = da1_dt1 + da1_dt2;
    %
%     a2 = (l1 / l2) * cos(t1 - t2);
%     da2_dt1 = (l1/l2) * (-sin(t1-t2)*w1); da2_dt2 = (l1/l2) * sin(t1-t2)*w2;
%     a2_dot = da2_dt1 + da2_dt2;
%     a2_dot = - l1/l2 * sin(t1 - t2) * (w1 -w2);
%    
da1_2_dt1_2 = - (l2/l1) * (m2/(m2+m1)) * (cos(t1-t2)*w1^2 + sin(t1-t2)*theta_s_ddot);
da1_2_dt2_2 = - (l2/l1) * (m2/(m2+m1)) * (cos(t1-t2)*w2^2 - sin(t1-t2)*theta_e_ddot);
da1_2_dt1_dt2 = (l2/l1) * (m2/(m2+m1)) * cos(t1-t2)*w1*w2;
a1_ddot = da1_2_dt1_2 + da1_2_dt2_2 + 2*da1_2_dt1_dt2;
%
da2_2_dt1_2 = - (l1 / l2) * (cos(t1-t2)*w1^2 + sin(t1-t2)*theta_s_ddot);
da2_2_dt2_2 = - (l1 / l2) * (cos(t1-t2)*w2^2 - sin(t1-t2)*theta_e_ddot);
da2_2_dt1_dt2 = (l1/l2) * cos(t1-t2)*w1*w2;
a2_ddot = da2_2_dt1_2 + da2_2_dt2_2 + 2*da2_2_dt1_dt2;
%
M_ddot = [0 a1_ddot
          a2_ddot 0];