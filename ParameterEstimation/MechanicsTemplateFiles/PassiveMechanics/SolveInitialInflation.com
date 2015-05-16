############ Set up output directory ##########################
read commands;SetOutput;

########### Read in the reference wall model #####################
read commands;ReadRefWallModel;
	
########### Read in the reference cavity model ###################
read commands;ReadRefCavityModel;
	
fem define mate;r;LV_CubicGuc
read commands;SolveInflation;

fem list elem;output_cavity_volume/LVCavityCurrent deformed total region $LV_CAVITY

############ Project registered surface data #######################
read commands;CalObjective;	
fem quit;
