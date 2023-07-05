function [ gamma_p, gamma_a, dgp_dln, dga_dln ] = comp_g_pa(l_n)
% Function to compute passive and length dependent gain
% and their respective partial derivative with respect to l_n
global k e_0m g Wp p

% % passive force gain
% if l_n <= 1,                gamma_p = 0;
%                             dgp_dln = 0;
% else,                       gamma_p = ( exp( k*( l_n - 1 ) / e_0m ) - 1 ) / ( exp(k) - 1 );
%                             dgp_dln = k*exp( k*( l_n - 1 ) / e_0m ) / ( e_0m*( exp(k) - 1 ) );
% end

Phip = zeros(p+1,1);
dPhip = zeros(p+1,1);

for j = 1:p+1,      Phip(j)  = l_n^(j-1);               end
for j = 2:p+1,      dPhip(j) = (j-1)*l_n^(j-2);         end

gamma_p = Wp*Phip;
dgp_dln = Wp*dPhip;

% active force gain
gamma_a = exp( -( ( l_n - 1 )^2 ) / g );
dga_dln = -2*gamma_a*( l_n - 1 ) / g;