function [ l_bar, dl_bar_dtheta, d2l_bar_d2theta  ] = comp_l_bar(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute muscle tendon complex length for all muscles

% works for both, geometry (2D) and geometry3 (3D)
global TRIlong TRIlat TRImed BIClong BICshort BRA S1 S2 S3 S4 S5 S6 urh S nm

l_bar           = NaN(nm, 1);
[ ~, d ] = size(TRIlong);
if d == 3,          J = comp_rot3(theta)';
else,               J = comp_rot(theta)';
end

TRIlong_g  = [ TRIlong(1:end-1,:);    (TRIlong(end,:)  - urh)*J + urh ]; 
TRIlat_g   = [ TRIlat(1:end-1,:);     (TRIlat(end,:)   - urh)*J + urh ];
TRImed_g   = [ TRImed(1:end-1,:);     (TRImed(end,:)   - urh)*J + urh ];
BIClong_g  = [ BIClong(1:end-1,:);    (BIClong(end,:)  - urh)*J + urh ];
BICshort_g = [ BICshort(1:end-1,:);   (BICshort(end,:) - urh)*J + urh ];
BRA_g      = [ BRA(1:end-1,:);        (BRA(end,:)      - urh)*J + urh ];


[ n, ~ ] = size(TRIlong_g); 
l_bar(1) = ones(1,n-1)*sqrt( diag( S1*TRIlong_g*TRIlong_g'*S1' ) );

[ n, ~ ] = size(TRIlat_g); 
l_bar(2) = ones(1,n-1)*sqrt( diag( S2*TRIlat_g*TRIlat_g'*S2' ) );

[ n, ~ ] = size(TRImed_g); 
l_bar(3) = ones(1,n-1)*sqrt( diag( S3*TRImed_g*TRImed_g'*S3' ) );

[ n, ~ ] = size(BIClong_g); 
l_bar(4) = ones(1,n-1)*sqrt( diag( S4*BIClong_g*BIClong_g'*S4' ) );

[ n, ~ ] = size(BICshort_g); 
l_bar(5) = ones(1,n-1)*sqrt( diag( S5*BICshort_g*BICshort_g'*S5' ) );

[ n, ~ ] = size(BRA_g); 
l_bar(6) = ones(1,n-1)*sqrt( diag( S6*BRA_g*BRA_g'*S6' ) );


%
z1         = TRIlong_g(end-1,:);            z2 = TRIlong_g(end,:);
dz2_dtheta = (TRIlong(end,:) - urh)*S*J;
dl_bar_dtheta(1) = -( z1 - z2 )*dz2_dtheta' /norm( z1 - z2 );
d2l_bar_d2theta(1) = dz2_dtheta*dz2_dtheta' /norm( z1 - z2 ) - ( z1 - z2 )*S'*dz2_dtheta' /norm( z1 - z2 ) - ( ( ( z1 - z2 )*dz2_dtheta' )^2 )/norm( z1 - z2 )^3;

%
z1 = TRIlat_g(end-1,:);             z2 = TRIlat_g(end,:);
dz2_dtheta = (TRIlat(end,:) - urh)*S*J;
dl_bar_dtheta(2) = -( z1 - z2 )*dz2_dtheta' /norm( z1 - z2 );
d2l_bar_d2theta(2) = dz2_dtheta*dz2_dtheta' /norm( z1 - z2 ) - ( z1 - z2 )*S'*dz2_dtheta' /norm( z1 - z2 ) - ( ( ( z1 - z2 )*dz2_dtheta' )^2 )/norm( z1 - z2 )^3;

%
z1 = TRImed_g(end-1,:);             z2 = TRImed_g(end,:);
dz2_dtheta = (TRImed(end,:) - urh)*S*J;
dl_bar_dtheta(3) = -( z1 - z2 )*dz2_dtheta' /norm( z1 - z2 );
d2l_bar_d2theta(3) = dz2_dtheta*dz2_dtheta' /norm( z1 - z2 ) - ( z1 - z2 )*S'*dz2_dtheta' /norm( z1 - z2 ) - ( ( ( z1 - z2 )*dz2_dtheta' )^2 )/norm( z1 - z2 )^3;
%
z1 = BIClong_g(end-1,:);            z2 = BIClong_g(end,:);
dz2_dtheta = (BIClong(end,:) - urh)*S*J;
dl_bar_dtheta(4) = -( z1 - z2 )*dz2_dtheta' /norm( z1 - z2 );
d2l_bar_d2theta(4) = dz2_dtheta*dz2_dtheta' /norm( z1 - z2 ) - ( z1 - z2 )*S'*dz2_dtheta' /norm( z1 - z2 ) - ( ( ( z1 - z2 )*dz2_dtheta' )^2 )/norm( z1 - z2 )^3;
%
z1 = BICshort_g(end-1,:);           z2 = BICshort_g(end,:);
dz2_dtheta = (BICshort(end,:) - urh)*S*J;
dl_bar_dtheta(5) = -( z1 - z2 )*dz2_dtheta' /norm( z1 - z2 );
d2l_bar_d2theta(5) = dz2_dtheta*dz2_dtheta' /norm( z1 - z2 ) - ( z1 - z2 )*S'*dz2_dtheta' /norm( z1 - z2 ) - ( ( ( z1 - z2 )*dz2_dtheta' )^2 )/norm( z1 - z2 )^3;
%
z1 = BRA_g(end-1,:);                z2 = BRA_g(end,:);
dz2_dtheta = (BRA(end,:) - urh)*S*J;
dl_bar_dtheta(6) = -( z1 - z2 )*dz2_dtheta' /norm( z1 - z2 );
d2l_bar_d2theta(6) = dz2_dtheta*dz2_dtheta' /norm( z1 - z2 ) - ( z1 - z2 )*S'*dz2_dtheta' /norm( z1 - z2 ) - ( ( ( z1 - z2 )*dz2_dtheta' )^2 )/norm( z1 - z2 )^3;
%

end
