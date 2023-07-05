
global TRIlong TRIlat TRImed BIClong BICshort BRA ...
       ground humerus urh hand S1 S2 S3 S4 S5 S6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Bodies in global coordinates
ground = [0 0];
humerus = [-0.017545 -0.007];
urh = [-0.011445 -0.2974];
hand = [-0.0011 -0.23559] + urh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Muscle Attachments - Initial Global: theta = 0 

TRIlong_P1 = [-0.0537 -0.0137] + ground;
TRIlong_P2 = [-0.0271 -0.1144] + humerus;
TRIlong_P3 = [-0.0318 -0.2264] + humerus;
TRIlong_P4 = [-0.0174 -0.2676] + humerus;
TRIlong_P5 = [-0.0219 0.0105] + urh;
TRIlong = [TRIlong_P1; TRIlong_P2; TRIlong_P3; TRIlong_P4; TRIlong_P5];

TRIlat_P1 = [-0.0060 -0.1265] + humerus;
TRIlat_P2 = [-0.0234 -0.1453] + humerus;
TRIlat_P3 = [-0.0318 -0.2264] + humerus;
TRIlat_P4 = [-0.0174 -0.2676] + humerus;
TRIlat_P5 = [-0.0219 0.0105] + urh;
TRIlat = [TRIlat_P1; TRIlat_P2; TRIlat_P3; TRIlat_P4; TRIlat_P5];

TRImed_P1 = [-0.0084 -0.1370] + humerus;
TRImed_P2 = [-0.0260 -0.1514] + humerus;
TRImed_P3 = [-0.0318 -0.2264] + humerus;
TRImed_P4 = [-0.0174 -0.2676] + humerus;
TRImed_P5 = [-0.0219 0.0105] + urh;
TRImed = [TRImed_P1; TRImed_P2; TRImed_P3; TRImed_P4; TRImed_P5];

BIClong_P1 = [-0.0392 0.0035] + ground;
BIClong_P2 = [-0.0289 0.0139] + ground;
BIClong_P3 = [0.0213 0.0179] + humerus;
BIClong_P4 = [0.0238 -0.0051] + humerus;
BIClong_P5 = [0.0135 -0.0283] + humerus;
BIClong_P6 = [0.0107 -0.0774] + humerus;
BIClong_P7 = [0.0170 -0.1213] + humerus;
BIClong_P8 = [0.0228 -0.1754] + humerus;
BIClong_P9 = [0.0075 -0.0484] + urh;
BIClong = [BIClong_P1; BIClong_P2; BIClong_P3; BIClong_P4; BIClong_P5; ...
           BIClong_P6; BIClong_P7; BIClong_P8; BIClong_P9];

BICshort_P1 = [0.0047 -0.0123] + ground;
BICshort_P2 = [-0.0071 -0.0400] + ground;
BICshort_P3 = [0.0112 -0.0758] + humerus;
BICshort_P4 = [0.0170 -0.1213] + humerus;
BICshort_P5 = [0.0228 -0.1754] + humerus;
BICshort_P6 = [0.0075 -0.0484] + urh;
BICshort = [BICshort_P1; BICshort_P2; BICshort_P3; BICshort_P4; ...
            BICshort_P5; BICshort_P6];

BRA_P1 = [0.0068 -0.1739] + humerus;
BRA_P2 = [-0.0032 -0.0239] + urh;
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