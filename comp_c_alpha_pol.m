function [c_alpha, dcos_dl] = comp_c_alpha_pol(l)

global Wc L_0 nm

p = 16; p = p - 1;
Phi = zeros(p+1,1); d_phi_dl = zeros(p+1,1); 



c_alpha = ones(nm, 1);     dcos_dl = zeros(nm, 1);

for i = 1:nm-3
    ln = l(i) / L_0(i);
    for j = 1:p+1,    Phi(j)        = ln^(j-1);         end
    for j = 2:p+1,    d_phi_dl(j)   = (j-1)*ln^(j-2);         end

    c_alpha(i) = Wc(i, :)*Phi;
    dcos_dl(i) = Wc(i, :)*d_phi_dl;
end







