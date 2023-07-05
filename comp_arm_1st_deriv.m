function g_dot = comp_arm_1st_deriv(theta, theta_dot, theta_ddot)
% Compute first order derivatives of skeletal dynamics - verified

global m1 m2 l1 l2 gr

t1 = theta(1);                  t2 = theta(2);
w1 = theta_dot(1);              w2 = theta_dot(2);
theta_s_ddot = theta_ddot(1);   theta_e_ddot = theta_ddot(2);    
%
C       = [0 (l2 / l1)*(m2/(m1+m2))*sin(t1-t2)
            -(l1/l2)*sin(t1-t2) 0];  
C_dot   = [0 (l2 / l1) * (m2 / (m2 + m1))*cos(t1-t2)*(w1 - w2)
            -(l1/l2)*cos(t1-t2)*(w1 - w2) 0];
        
%     gt      = [gr/l1*sin(t1) gr/l2*sin(t2)]';                               
gt_dot  = [(gr/l1)*cos(t1)*w1 (gr/l2)*cos(t2)*w2]';
%
g_dot = C_dot * [w1^2 w2^2]' + C * [2*w1*theta_s_ddot 2*w2*theta_e_ddot]' + gt_dot ;

    
end