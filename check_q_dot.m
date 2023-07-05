function [q_dot] = check_q_dot(q_dot, q)
% Constrain activation dynamics between q_min (0.01) and q_max (1)
global nm q_min q_max

for i = 1:nm
   if (q(i)>=q_max)&&(q_dot(i)>0)
       q_dot(i) = 0;
   elseif (q(i)<=q_min)&&(q_dot(i)<0)
       q_dot(i) = 0;
   end  
end

end