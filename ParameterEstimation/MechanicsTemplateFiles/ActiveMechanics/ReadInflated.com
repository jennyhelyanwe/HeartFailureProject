#################### Read in the model at the end of inflation ################
#$NAME="LV_Inflation_OptC1";
#$FILE=$output.${NAME};
$FILE="WarmStartSolution/CurrentContracted_1";
fem define initial;r;$FILE region $WALL;
fem define mapping;r;LV_CubicMapAll region $WALL;
fem define solve;r;LV_Cubic region $WALL;

fem evaluate residual;EndInf_Opt wrt geom_params region $WALL;
system('tail EndInf_Opt.opresi');

########################### Compute LV Volume #################################
read commands;Extract_Pressure;
read commands;Cal_CavityVolumeMF;

#############################################################################################################

########## Redefine files required for the wall mechanics simulation ############################
fem define initial;r;$FILE region $WALL;
fem define mapping;r;LV_CubicMapAll region $WALL;
fem define solve;r;LV_Cubic region $WALL;
fem define equation;r;LV_Cubic region $WALL;

fem export node;LVEndOfInflation_OptC1 as LVEndOfInflation_OptC1 field region $WALL;
fem export elem;LVEndOfInflation_OptC1 as LVEndOfInflation_OptC1 field region $WALL;

$NAME="LVModel_1";
$FILE=$outputAll.${NAME};
fem export node;$FILE as field LVModel region $WALL;
fem export elem;$FILE as field LVModel region $WALL;
