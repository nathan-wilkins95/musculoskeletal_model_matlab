function [y] = tanh_int(x)

y   =      (1825185783670951*log(exp(75*x)*exp(-81/80) + 1))/13194139533312000 + ...
            (1825185783670951*log(exp(150*x)*exp(-81/40) - exp(75*x)*exp(-81/80) + exp(300*x)*exp(-81/20) - ...
            exp(225*x)*exp(-243/80) + 1))/13194139533312000;

end