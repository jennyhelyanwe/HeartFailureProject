	

	read commands;SetOutput;
	read commands;ReadRefWallModel;
	read commands;ReadRefCavityModel;


	#################### Read in the model at the end of inflation ################
    	$NAME="LV_Inflation_OptC1";
    	$FILE=$output.${NAME};
	fem define initial;r;$FILE region $WALL;
	fem define mapping;r;LV_CubicMapAll region $WALL;  
	fem define solve;r;LV_Cubic region $WALL;

	fem evaluate residual;EndInf_Opt wrt geom_params region $WALL;
	system('tail EndInf_Opt.opresi');
	
	########################### Compute LV Volume #################################
	$ca_file="TCa_ZERO"
	$LINE_TCa=16;
	$TCa=`sed -e s%D%E% $ca_file."ipacti" | awk -v line=$LINE_TCa 'NR==line{printf("%.2f",\$9)}'` # Extract Calcium value
	read commands;Extract_Pressure;
	read commands;Cal_CavityVolumeMF;

	MySubs::writeCalciumSystole($TCa,$LV_P,$LV_CAVITY_VOLUME);
	
	#############################################################################################################

	########## Redefine files required for the wall mechanics simulation ########
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

	fem quit;
