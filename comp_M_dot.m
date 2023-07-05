function [M_dot, M_dot_inv] = comp_M_dot(theta, theta_dot)
% compute time derivative of inertia matrix M and the corresponding inverse
% M_inv
global m1 m2 l1 l2

t1 = theta(1);  t2 = theta(2);  w1 = theta_dot(1);  w2 = theta_dot(2);

M_dot = [0  -(l2 / l1) * (m2 / (m1 + m2))*sin(t1 - t2)*(w1-w2)
         -l1/l2*sin(t1 - t2)*(w1 -w2) 0];
%
a1 = (l2 / l1) * (m2 / (m1 + m2)) * cos(t1 - t2);
da1_dt1 = - (l2 / l1) * (m2 / (m1 + m2)) * sin(t1 - t2) * w1; da1_dt2 = (l2 / l1) * (m2 / (m1 + m2)) * sin(t1 - t2) * w2;
a1_dot = da1_dt1 + da1_dt2;
%
a2 = (l1 / l2) * cos(t1 - t2);
da2_dt1 = (l1/l2) * (-sin(t1-t2)*w1); da2_dt2 = (l1/l2) * sin(t1-t2)*w2;
a2_dot = da2_dt1 + da2_dt2;
%
M_dot_inv = 1/((1-a1*a2)^2) * [a1_dot*a2+a1*a2_dot  -a1_dot*(1-a1*a2)-a1*(a1_dot*a2+a1*a2_dot)
                                -a2_dot*(1-a1*a2)-a2*(a1_dot*a2+a1*a2_dot)     a1_dot*a2+a1*a2_dot];