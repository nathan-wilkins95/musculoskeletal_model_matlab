function [r] = sigmoid_c(x)

a = 10; b = 0.5;

r = 1 ./ (1 + exp(-a*(x - b)));

end