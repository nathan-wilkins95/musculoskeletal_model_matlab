function [theta_des_dot, theta_des_ddot, theta_des_3dot, theta_des_4dot] = comp_des_bell(t)

global theta_sf theta_si theta_ef theta_ei T d2r

theta_dfi = [(theta_sf - theta_si) (theta_ef - theta_ei)]';

if (t > T) && (t < 2*(T+0.01))

    theta_des_dot       = -( (-60*(t-T)^3)/T^4  + 30*(t-T)^4/T^5 + 30*(t-T)^2/T^3 )*theta_dfi*d2r;
    theta_des_ddot      = -( (-180*(t-T)^2)/T^4  + 120*(t-T)^3/T^5 + 60*(t-T)/T^3 )*theta_dfi*d2r;
    theta_des_3dot      = -( (-360*(t-T))/T^4  + 360*(t-T)^2/T^5 + 60/T^3 )*theta_dfi*d2r;
    theta_des_4dot      = -( (-360)/T^4  + 720*(t-T)/T^5 )*theta_dfi*d2r;

elseif t <= T
    
    theta_des_dot       = ( (-60*t^3)/T^4  + 30*t^4/T^5 + 30*t^2/T^3 )*theta_dfi*d2r;
    theta_des_ddot      = ( (-180*t^2)/T^4  + 120*t^3/T^5 + 60*t/T^3 )*theta_dfi*d2r;
    theta_des_3dot      = ( (-360*t)/T^4  + 360*t^2/T^5 + 60/T^3 )*theta_dfi*d2r;
    theta_des_4dot      = ( (-360)/T^4  + 720*t/T^5 )*theta_dfi*d2r;

end
    

end