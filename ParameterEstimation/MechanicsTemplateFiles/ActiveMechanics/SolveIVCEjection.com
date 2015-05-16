

##### Auxiliary parameters
$ITERS=10;
$ERROR_TOLERANCE=1e-6;
$Pee=19.0;

## Get the preload

read commands;Extract_Pressure;
################# Start with Phase 2: IVC ####################################
$No_Steps=30;
$PIncr=($Pee-$LV_P)/$No_Steps;
$dis_step=1.0/$No_Steps;

##################### Set up calcium level ##############################
# This perl script contains a series of calcium concentration at given time, 
# the entire simulation will be driven by the calcium levels and will stop
# when the last calcium level is read
$index=1;
$output_index=1;
$continueReading=1; # continueReading = false if the last calcium reading was reached
$TCalVIndex=1; # Calculate the number of times cavity volume was calculated. 
$Ca=0;



		system('perl CopyInitial.pl');
		set echo;
		#==================================================================================================
		print "\033[0;30;42m ============================================================ \033[0m\n"; 
		print "\033[0;30;42m                    Solution Method                           \033[0m\n"; 
		print "\033[0;30;42m ============================================================ \033[0m\n"; 
		#==================================================================================================
		set echo;
		fem define initial;r;LV_CubicPreEpiBase region $WALL;
		fem define mapping;r;LV_CubicMapAll region $WALL;
		fem define solve;r;LV_Cubic region $WALL;  
		
	
while ($LV_P<$Pee)
{

		########################################################################		
		############################ Read in Calcium ###########################
		########################################################################
	
	      	
	    ################ Displacement Boundary Conditions #####################
		# Perform displacement 
		# Solve with defined increment to apply disired boundary conditions
		print "\033[0;30;43m    Performing Displacement Boundary Conditions       \033[0m\n";
      		fem solve increment $dis_step iter $ITERS  error $ERROR_TOLERANCE;
      		fem update pressure boundary elements WALL auxillary 1 increment -$dis_step;
      		fem solve increment 0.0 iter $ITERS error $ERROR_TOLERANCE region $WALL;
		read commands;Extract_Pressure;
		########################################################################

		####################Incremental Approach #######################
 		# Create ipacti file with the new calcium levle which is constant spatially
		#if ($index == 1)
		#{		
		#	$TCa=0;
		#} else
 		#{
 			#$Ca=$Ca+0.1;			
			#$TCa=$Ca**3/($Ca**3+0.5**3)*100;
			$TCa=$TCa+2;
 		#}
		MySubs::changeCalciumEspen_TCa($TCa,$index);
		$TCa_file="calcium0xx_".${index};
		fem define active;r;$TCa_file region $WALL;
		print "\033[0;30;43m     Increase TCa Level to  ${TCa}      \033[0m\n"

		fem solve increment 0.0 iter $ITERS error $ERROR_TOLERANCE region $WALL;
	    	read commands;Extract_Pressure;
		########################################################################
		
		#######################################################################
		############################ Increase Pressure ########################
		#######################################################################
		$LV_P_Current=$LV_P+$PIncr;
		
		print "\033[0;30;43m     Increase Pressure Level to  ${LV_P_Current}      \033[0m\n"
		# Update the pressure and correct the pressure
		fem update pressure boundary elements WALL auxillary 1 increment $PIncr	
		# Solve with 0 pressure increment to ensure system is converged ::: Maybe optional
		fem solve increment 0.0 iter $ITERS error $ERROR_TOLERANCE region $WALL;
		read commands;Extract_Pressure;
		
		######## Write out current state 
		fem def init;w;currentInitial region $WALL;
	
		######## Compute the cavity volume
		read commands;Cal_CavityVolume;
		fem def init;r;currentInitial region $WALL;
		fem def solve;r;LV_Cubic region $WALL;
		fem def mapping;r;LV_CubicMapAll region $WALL;
		

		$NAME="LV_Systole_FastTrack_".${output_index};
		$FILE=$output.${NAME};
		MySubs::writeCalciumSystole($TCa,$LV_P,$index,$LV_CAVITY_VOLUME);
		$FILE2="outputAll/LVModel_".${output_index};
		fem export nodes;$FILE2 field as LVModel region $WALL;
		fem export eleme;$FILE2 field as LVModel region $WALL;
		########################## Output ###############################
		# Save the current status as the initial condition for the next phase
		fem define initial;w;$FILE region $WALL;
		
		# Export the model first
		fem export nodes;$FILE field as LVSystole_FastTrack region $WALL;
		fem export elements;$FILE field as LVSystole_FastTrack region $WALL;
			
		fem evaluate residual;TestingIVC wrt geom_params region $WALL;
				
		# Update extension ratio and strain
		fem update gauss strain extension_ratios region $WALL;
		fem export gauss;${FILE}."_gauss_ER_FastTrack" yg as gauss_strain region $WALL;
		fem update gauss strain region $WALL;
		fem export gauss;${FILE}."_gauss_strain_FastTrack" yg as gauss_strain region $WALL;
		fem update gauss strain wall region $WALL;
		fem export gauss;${FILE}."_gauss_wall_strain_FastTrack" yg as gauss_wall_strain region $WALL;
		fem update gauss stress total cauchy region $WALL;
		fem export gauss;${FILE}."_gauss_stress_FastTrack" yg as gauss_stress region $WALL;
		fem update gauss stress passive cauchy region $WALL;
		fem export gauss;${FILE}."_passive_gauss_stress_FastTrack" yg as gauss_stress region $WALL;
		fem update gauss stress active cauchy region $WALL;
		fem export gauss;${FILE}."_active_gauss_stress_FastTrack" yg as gauss_stress region $WALL;
		#########################################################################
		$index=$index+1;
		$output_index=$output_index+1;
		#==========================================================================================================
		print "\033[0;30;42m ============================================================ \033[0m\n"; 
		print "\033[0;30;42m            Going to Next Ca/Pressure Activation Level        \033[0m\n"; 
		print "\033[0;30;42m ============================================================ \033[0m\n"; 
		#==========================================================================================================
	
}

## Write out the Tca at ES;
MySubs::changeCalciumEspen_ESTCa($TCa);

	      		
	
	
