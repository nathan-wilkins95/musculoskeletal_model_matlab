
global TRIlong TRIlat TRImed BIClong BICshort BRA ...
       ground humerus urh hand S1 S2 S3 S4 S5 S6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Bodies in global coordinates
ground = [0 0 0];
humerus = [-0.017545 -0.007 0.17];  clavicle = [-0.011096 0.0063723 0.054168];
urh = [-0.011445 -0.2974 0.154635]; scapula = [-0.054694 -0.035032 -0.03734];
hand = [-0.0011 -0.23559 0.1577] + urh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Muscle Attachments - Initial Global: theta = 0 

DELT1_P1 = [0.00896 -0.11883 0.00585]   + humerus;
DELT1_P2 = [0.01623 -0.11033 0.00412]   + humerus;
DELT1_P3 = [0.04347 -0.03252 0.00099]   + clavicle;
DELT1_P4 = [-0.014 0.01106 0.08021]     + scapula;

TRIlong_P1 = [-0.05365 -0.01373 0.1473] + ground;
TRIlong_P2 = [-0.0271 -0.1144 -0.00664] + humerus;
TRIlong_P3 = [-0.0318 -0.2264 -0.01217] + humerus;
TRIlong_P4 = [-0.0174 -0.2676 -0.01208] + humerus;
TRIlong_P5 = [-0.0219 0.0105 -0.00078] + urh;
TRIlong = [TRIlong_P1; TRIlong_P2; TRIlong_P3; TRIlong_P4; TRIlong_P5];

TRIlat_P1 = [-0.0060 -0.1265 0.00428] + humerus;
TRIlat_P2 = [-0.0234 -0.1453 0.00928] + humerus;
TRIlat_P3 = [-0.0318 -0.2264 -0.01217] + humerus;
TRIlat_P4 = [-0.0174 -0.2676 -0.01208] + humerus;
TRIlat_P5 = [-0.0219 0.0105 -0.00078] + urh;
TRIlat = [TRIlat_P1; TRIlat_P2; TRIlat_P3; TRIlat_P4; TRIlat_P5];

TRImed_P1 = [-0.0084 -0.1370 -0.00906] + humerus;
TRImed_P2 = [-0.0260 -0.1514 -0.0108] + humerus;
TRImed_P3 = [-0.0318 -0.2264 -0.01217] + humerus;
TRImed_P4 = [-0.0174 -0.2676 -0.01208] + humerus;
TRImed_P5 = [-0.0219 0.0105 -0.00078] + urh;
TRImed = [TRImed_P1; TRImed_P2; TRImed_P3; TRImed_P4; TRImed_P5];

BIClong_P1 = [-0.039235 0.00347 0.14795] + ground;
BIClong_P2 = [-0.028945 0.01391 0.15639] + ground;
BIClong_P3 = [0.02131 0.01793 0.01028] + humerus;
BIClong_P4 = [0.02378 -0.00511 0.01201] + humerus;
BIClong_P5 = [0.01345 -0.02827 0.00136] + humerus;
BIClong_P6 = [0.01068 -0.07736 -0.00165] + humerus;
BIClong_P7 = [0.01703 -0.12125 0.00024] + humerus;
BIClong_P8 = [0.0228 -0.1754 -0.0063] + humerus;
BIClong_P9 = [0.00751 -0.04839 0.02179] + urh;
BIClong = [BIClong_P1; BIClong_P2; BIClong_P3; BIClong_P4; BIClong_P5; ...
           BIClong_P6; BIClong_P7; BIClong_P8; BIClong_P9];

BICshort_P1 = [0.0047 -0.0123 0.13475] + ground;
BICshort_P2 = [-0.0071 -0.0400 0.14507] + ground;
BICshort_P3 = [0.0112 -0.0758 -0.01101] + humerus;
BICshort_P4 = [0.0170 -0.1213 -0.01079] + humerus;
BICshort_P5 = [0.0228 -0.1754 -0.0063] + humerus;
BICshort_P6 = [0.0075 -0.0484 0.02179] + urh;
BICshort = [BICshort_P1; BICshort_P2; BICshort_P3; BICshort_P4; ...
            BICshort_P5; BICshort_P6];

BRA_P1 = [0.0068 -0.1739 -0.0036] + humerus;
BRA_P2 = [-0.0032 -0.0239 0.0009] + urh;
BRA = [BRA_P1; BRA_P2];
                    
%%
S1 = comp_myS(TRIlong);         S2 = comp_myS(TRIlat);
S3 = comp_myS(TRImed);          S4 = comp_myS(BIClong);
S5 = comp_myS(BICshort);        S6 = comp_myS(BRA);

%% Delete superfluous parameters
clear TRIlong_P1 TRIlong_P2 TRIlong_P3 TRIlong_P4 TRIlong_P5
clear TRIlat_P1 TRIlat_P2 TRIlat_P3 TRIlat_P4 TRIlat_P5
clear TRImed_P1 TRImed_P2 TRIlong_P3 TRImed_P4 TRImed_P5
clear BIClong_P1 BIClong_P2 BIClong_P3 BIClong_P4 BIClong_P5 BIClong_P6 BIClong_P7 BIClong_P8 BIClong_P9
clear BICshort_P1 BICshort_P2 BICshort_P3 BICshort_P4 BICshort_P5 BICshort_P6
clear BRA_P1 BRA_P2