function [tau_fe] = tau_to_tau_fe(tau)

c = comp_c(tau);

tau_f  = subplus( tau) + c;
tau_e  = subplus(-tau) + c;

tau_fe = [ tau_f(1); -tau_e(1); tau_f(2); -tau_e(2) ];