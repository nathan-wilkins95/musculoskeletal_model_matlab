function [y] = tanh_fun(x)

a = 300;        c = 51.8749;         d = 0.0135;     as = 0.25;     
y = c*( exp(a*(x - d)) )./( exp(a*(x - d)) + exp(-as*a*(x - d)) );

end