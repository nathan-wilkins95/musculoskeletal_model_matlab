function [ ln_dot, dh, flg ] = ncomp_ln_dot( gamma, q )

global p1 p2 p3 p4 p5 p6 p7

ga = gamma(1);  gc = gamma(2);

if( ( gc > 0 ) && ( gc < q*ga*p3*p4 ) )
    %
    if gc <= q*ga
        ln_dot  = ( p1 + p2*q )*( gc - q*ga ) / ( q*ga + gc/p5 );                flg = 1;
        dl_dq   = p2*( gc - q*ga ) / ( q*ga + gc/p5 ) - ( p1 + p2*q )*ga / ( q*ga + gc/p5 ) - ( p1 + p2*q )*( gc - q*ga )*ga / ( q*ga + gc/p5 )^2;
        dl_dga  = -q*( p1 + p2*q )*(1+1/p5)*gc  / ( q*ga + gc/p5 )^2;
        dl_dgc  =  q*( p1 + p2*q )*(1+1/p5)*ga  / ( q*ga + gc/p5 )^2;
    else
        ln_dot  = ( p1 + p2*q )*( gc - q*ga )*p6 / ( p7*( q*ga*p3 - gc ) );      flg = 2;
        dl_dq   = p2*( gc - q*ga )*p6 / ( p7*( q*ga*p3 - gc ) ) - ( p1 + p2*q )*ga*p6 / ( p7*( q*ga*p3 - gc ) )...
                  - ( p1 + p2*q )*( gc - q*ga )*p6*ga*p3 / ( p7*( q*ga*p3 - gc )^2 );
        dl_dga  = -q*( p1 + p2*q )*( p6 / p7 )*gc*( p3 - 1 )  / ( q*ga*p3 - gc )^2;
        dl_dgc  =  q*( p1 + p2*q )*( p6 / p7 )*ga*( p3 - 1 )  / ( q*ga*p3 - gc )^2;
    end
    %
else
   %
   if gc <= 0
        ln_dot = ( p1 + p2*q )*( (1+1/p5)*gc - q*ga )/( q*ga );                      flg  = 3;
        dl_dq  = p2*( (1+1/p5)*gc - q*ga )/( q*ga ) - ( p1/q + p2 ) - ( p1 + p2*q )*( (1+1/p5)*gc - q*ga )/( ga*q^2 );
        dl_dga = -( p1 + p2*q )*( 1 + 1/p5 )*gc / ( q*ga^2 );
        dl_dgc =  ( p1 + p2*q )*( 1 + 1/p5 )*ga / ( q*ga^2 );
   else
        c1 = p1*p6/( p7*p3*(1-p4) );                           c2 = p2*p6/( p7*p3*(1-p4) );
        c3 = (p3*p4 - 1)*( 1 - 2*p4)/( 1 - p4);                c4 = p3*p4;
        c5 = (p3*p4 - 1 )/(p3*( 1 - p4) );

        ln_dot = ( c1 + c2*q )*( c3 - gc + c4*q*ga + c5*gc/(q*ga) );          flg  = 4;
        dl_dq  = c2*( c3 - gc + c4*q*ga + c5*gc/(q*ga) ) + ( c1 + c2*q )*( c4*ga - c5*gc/(ga*q^2) );
        dl_dga = ( c1 + c2*q )*( c4*q - c5*gc/( q*ga^2) );
        dl_dgc = ( c1 + c2*q )*( -1 + c5/(q*ga) );
   end
end

dh = [ dl_dq dl_dga dl_dgc ]';

end
