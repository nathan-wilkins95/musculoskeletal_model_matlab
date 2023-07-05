function [Xi_p, Xi_p_dot] = comp_Xi_p_dot(Xi, Xi_dot) 

global nm

Xi_p        = pinv(Xi); 

Xi_p_dot    = -Xi_p*Xi_dot*Xi_p + Xi_p*Xi_p'*Xi_dot'*(eye(nm) - Xi*Xi_p) + (eye(nm) - Xi_p*Xi)*Xi_dot'*(Xi_p'*Xi_p);