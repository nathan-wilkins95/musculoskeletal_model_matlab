function [y] = sigmoid(x)

a = 500;     c = 51.8749;         d = 0.016;
y = c./(1 + exp(-a*(x - d) ));

end