function [norm_rate] = Prochazka_Ia_rates(l, lm_dot)
%     """ Compute Prochazka Ia rates """
global nm L_0 p8
a=4.3; b=2; c=10;

norm_rate = zeros(nm, 1);

for i = 1:nm
    opt_l = L_0(i) * 1000;
    max_v = p8 * opt_l;
    fiber_l = l(i) * 1000;
    fiber_v = lm_dot(i) * 1000;
    rate = a * sign(fiber_v) * exp(0.6 * log(max(min(abs(fiber_v), max_v), 0.01))) + b * (min(fiber_l, 1.5 * opt_l) - opt_l) + c;
    norm_rate(i) = max(rate / (a * exp(0.6 * log(max_v)) + b * 0.5 * opt_l + c), 0);
end
end