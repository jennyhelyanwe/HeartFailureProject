	#################### Read in the model at the end of inflation ################
    	
	#fem define initial;r;${OUTPUT_DIR}."EndOfInfla" region $WALL
	fem define initial;r;InitialInflated region $WALL;
	fem define mapping;r;LV_CubicMapAll region $WALL;  
	fem define solve;r;LV_Cubic region $WALL;


	########################### Compute LV Volume #################################
	read commands;Extract_Pressure;
	read commands;Cal_CavityVolume;


	#############################################################################################################

	########## Redefine files required for the wall mechanics simulation ############################
	fem define initial;r;InitialInflated region $WALL;
	fem define mapping;r;LV_CubicMapAll region $WALL;  
	fem define solve;r;LV_Cubic region $WALL;
	fem define equation;r;LV_Cubic region $WALL;
	
	fem export node;LVEndOfInflation as LVEndOfInflation field region $WALL;
	fem export elem;LVEndOfInflation as LVEndOfInflation field region $WALL;

	## Calculate deformed coordinates of the surface points
	fem define data;r;BackTransformedUPFSurfaceRotated_ModelTree_Epi_ED region $WALL;	
	fem group faces 5,9,13,16,21,25,29,32,37,41,45,48,53,57,61,64,51,55,59,62 as EPI region $WALL;
	fem def xi;c closest_face faces EPI region $WALL;
	
	## Calculate the deformed coordinates
	fem define data;c from_xi deformed region $WALL;
	fem define data;w;UPFEpiED_ModelPred region $WALL;
	fem export data;UPFEpiED_ModelPred as UPFEpiED_ModelPred region $WALL;
	
	
	## Calculate deformed coordinates of the surface points
	fem define data;r;BackTransformedUPFSurfaceRotated_ModelTree_Endo_ED region $WALL;	
	fem group faces 4,8,12,15,20,24,28,31,36,40,44,47,52,56,60,63,51,55,59,62 as ENDO region $WALL;
	fem def xi;c closest_face faces ENDO region $WALL;	
	
	## Calculate the deformed coordinates
	fem define data;c from_xi deformed region $WALL;
	fem define data;w;UPFEndoED_ModelPred region $WALL;
	fem export data;UPFEndoED_ModelPred as UPFEndoED_ModelPred region $WALL;
	
	
