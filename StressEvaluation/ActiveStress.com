# This script reads in a model solution (at maximum TCa) and evaluates the principle
# stress at the gauss points. 
fem def para;r;LV_Cubic
fem def coor;r;LV_Cubic
fem def base;r;LV_Cubic

fem def node;r;DSWall
fem def elem;r;DSWall

fem def fibre;r;DTMRI_CIMModel
fem def elem;r;DTMRI_CIMModel fibre

fem def equa;r;LV_Cubic
fem def mate;r;LV_CubicGuc

fem def acti;r;ActivationOpt/TCa_12
fem def init;r;PostContraction_BC/CurrentContracted_12
fem def solve;r;LV_Cubic
fem list acti
fem solve increm 0.0 error 1e-6 iter 10
fem list stress;ActiveStress_12 active full

fem quit

