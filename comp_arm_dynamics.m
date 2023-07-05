function [w, g] = comp_arm_dynamics(theta, theta_dot)
% compute skeletal dynamics - verified
global m1 m2 l1 l2 gr

t1 = theta(1);      t2 = theta(2);
w1 = theta_dot(1);  w2 = theta_dot(2);    
    
C = [0 (l2 / l1)*(m2/(m1+m2))*sin(t1-t2)
    -(l1/l2)*sin(t1-t2) 0];  
gt = [gr/l1*sin(t1) gr/l2*sin(t2)]';                               

g = C * [w1^2 w2^2]' + gt;                                  
w = [w1 w2]';
  
end