
	
	############ Set up output directory ##########################
	read commands;SetOutput;

	########### Read in the reference wall model #####################
	read commands;ReadRefWallModel;
	
	########### Read in the reference cavity model ###################
	read commands;ReadRefCavityModel;
	      
    	read commands;ReadInflated;
	read commands;SolveIVCEjection;
	
	fem define initial;w;InitialEndSystole region $WALL;

	################# Output ######################
	$NAME="LV_EndSystole_Initial";
	$FILE=$output.${NAME};
	fem define initial;w;$FILE region $WALL;

	# Export the model first
	fem export nodes;$FILE field as LVEndSystole_Initial region $WALL;
	fem export elements;$FILE field as LVEndSystole_Initial region $WALL;
	## Update the fibre field
	fem update gauss deformed_fibres region $WALL;
	fem export gauss;${FILE}."_gauss_Fibre" yg as gauss_fibre region $WALL;

	######## Save the strains and stresses to output folder ################
    	fem update gauss strain extension_ratios region $WALL;
    	fem export gauss;${FILE}."_gauss_ER" yg as gauss_strain region $WALL;
    	fem update gauss strain region $WALL;
    	fem export gauss;${FILE}."_gauss_strain" yg as gauss_strain region $WALL;
    	fem update gauss stress total cauchy region $WALL;
   	fem export gauss;${FILE}."_gauss_stress" yg as gauss_stress region $WALL;
    	fem update gauss stress passive cauchy region $WALL;
    	fem export gauss;${FILE}."_passive_gauss_stress" yg as gauss_stress region $WALL;
    	fem update gauss stress active cauchy region $WALL;
    	fem export gauss;${FILE}."_active_gauss_stress" yg as gauss_stress region $WALL;
	
	read commands;CalObjective;
	fem quit;
