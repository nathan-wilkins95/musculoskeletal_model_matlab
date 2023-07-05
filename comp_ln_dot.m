function [ ln_dot, mflag, ln_dot2, fl, gl, mpartials ] = comp_ln_dot(gamma_c,gamma_a,q)
% compute normalized muscle velocity
global p1 p2 p3 p4 p5 p6 p7

if( ( gamma_c > 0 ) && ( gamma_c < q*gamma_a*p3*p4 ) )
    %
    if gamma_c <= q*gamma_a
        ln_dot  = ( p1 + p2*q )*( gamma_c - q*gamma_a ) / ( q*gamma_a + gamma_c/p5 );                mflag = 1;
        dl_dq   = p2*( gamma_c - q*gamma_a ) / ( q*gamma_a + gamma_c/p5 ) - ( p1 + p2*q )*gamma_a / ( q*gamma_a + gamma_c/p5 ) - ( p1 + p2*q )*( gamma_c - q*gamma_a )*gamma_a / ( q*gamma_a + gamma_c/p5 )^2;
        dl_dga  = -q*( p1 + p2*q )*(1+1/p5)*gamma_c  / ( q*gamma_a + gamma_c/p5 )^2;
        dl_dgc  =  q*( p1 + p2*q )*(1+1/p5)*gamma_a  / ( q*gamma_a + gamma_c/p5 )^2;
    else
        ln_dot  = ( p1 + p2*q )*( gamma_c - q*gamma_a )*p6 / ( p7*( q*gamma_a*p3 - gamma_c ) );      mflag = 2;
        dl_dq   = p2*( gamma_c - q*gamma_a )*p6 / ( p7*( q*gamma_a*p3 - gamma_c ) ) - ( p1 + p2*q )*gamma_a*p6 / ( p7*( q*gamma_a*p3 - gamma_c ) )...
                  - ( p1 + p2*q )*( gamma_c - q*gamma_a )*p6*gamma_a*p3 / ( p7*( q*gamma_a*p3 - gamma_c )^2 );
        dl_dga  = -q*( p1 + p2*q )*( p6 / p7 )*gamma_c*( p3 - 1 )  / ( q*gamma_a*p3 - gamma_c )^2;
        dl_dgc  =  q*( p1 + p2*q )*( p6 / p7 )*gamma_a*( p3 - 1 )  / ( q*gamma_a*p3 - gamma_c )^2;
    end
    %
else
   %
   if gamma_c <= 0
%        ( p1 + p2*q )
%        ( (1+1/p5)*gamma_c - q*gamma_a )
%        ( q*gamma_a )
        ln_dot = ( p1 + p2*q )*( (1+1/p5)*gamma_c - q*gamma_a )/( q*gamma_a );                      mflag  = 3;
        dl_dq  = p2*( (1+1/p5)*gamma_c - q*gamma_a )/( q*gamma_a ) - ( p1/q + p2 ) - ( p1 + p2*q )*( (1+1/p5)*gamma_c - q*gamma_a )/( gamma_a*q^2 );
        dl_dga  = -( p1 + p2*q )*( 1 + 1/p5 )*gamma_c / ( q*gamma_a^2 );
        dl_dgc  =  ( p1 + p2*q )*( 1 + 1/p5 )*gamma_a / ( q*gamma_a^2 );
   else
        c1 = p1*p6/( p7*p3*(1-p4) );                           c2 = p2*p6/( p7*p3*(1-p4) );
        c3 = (p3*p4 - 1)*( 1 - 2*p4)/( 1 - p4);                c4 = p3*p4;
        c5 = (p3*p4 - 1 )/(p3*( 1 - p4) );

        ln_dot = ( c1 + c2*q )*( c3 - gamma_c + c4*q*gamma_a + c5*gamma_c/(q*gamma_a) );          mflag  = 4;
        dl_dq  = c2*( c3 - gamma_c + c4*q*gamma_a + c5*gamma_c/(q*gamma_a) ) + ( c1 + c2*q )*( c4*gamma_a - c5*gamma_c/(gamma_a*q^2) );
        dl_dga  = ( c1 + c2*q )*( c4*q - c5*gamma_c/( q*gamma_a^2) );
        dl_dgc  = ( c1 + c2*q )*( -1 + c5/(q*gamma_a) );

        
   end
end

fl      = ln_dot - dl_dq*q;             gl = dl_dq;
ln_dot2 = fl + gl*q;  

ln_dot    = ln_dot2;
mpartials = [ dl_dq dl_dga dl_dgc ]';
end
