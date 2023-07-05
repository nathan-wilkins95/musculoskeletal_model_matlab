function [ x, x_dot, x_ddot ] = comp_sin( t, par )

a = par(1);          w = par(2);      
p = par(3);          o = par(4);      

x      = diag(a)*sin( w*t + p ) + o;
x_dot  = diag(a)*diag(w)*cos( w*t + p );
x_ddot = -diag(a)*diag(w)*diag(w)*sin( w*t + p );