function [M, M_inv] = comp_M(theta)
% compute inertia matrix M and its inverse M_inv
global m1 m2 l1 l2

a1 = (l2 / l1) * (m2 / (m1 + m2)) * cos(theta(1) - theta(2));
a2 = (l1 / l2) * cos(theta(1) - theta(2));

M = [1 a1
    a2 1];

M_inv = 1/(1-a1*a2)*[1 -a1
                     -a2 1];