########################################################################################
################################# Optimization #########################################
########################################################################################
use POSIX
# This program is designed to optimize the stiffness parameter C1 of the transversely-isotropic
# law using a set of 3D material points obtained from tag tracking. The target data points will
# be material points tracked at frame 6 (hence end-diastolic state), and the landmark points 
# are points from frame 48. 

########################### Set up initial files ######################################


	use POSIX;
	set echo on;
	

	############ Set up output directory ##########################
	read commands;SetOutput;

	########### Read in the reference wall model #####################
	read commands;ReadRefWallModel;
	
	########### Read in the reference cavity model ###################
	read commands;ReadRefCavityModel;
	
	
	#################### Read in pre-solved ED model ################
	
	#################### Solve iteratively ################
	$TCa_goal=`sed -e s%D%E% TCa_OptTCa_current.ipacti | awk -v line=16 'NR==line{printf("%.5f",\$7)}'`;
	$TCa_previous = `sed -e s%D%E% TCa_previous.ipacti | awk -v line=16 'NR==line{printf("%.5f",\$7)}'`;
	$TCa_increm = $TCa_goal - $TCa_previous;

	print "\n \033[0;30;42m ================================================================================== \033[0m\n";
	print "\033[0;30;42m       Solving iteratively from TCa = $TCa_previous to reach TCa = $TCa_goal    	\033[0m\n"; 
	print "\033[0;30;42m ================================================================================== \033[0m\n";

	
	#if (abs($TCa_increm)>20)
	#{
	#		$MAXIMUM_INCREM = 10;
	#} elsif (abs($TCa_increm) > 10){
	#		$MAXIMUM_INCREM = 8;
	#} elsif (abs($TCa_increm) > 5){
	#		$MAXIMUM_INCREM = 5; 
	#} elsif (abs($TCa_increm) > 3) {
	#		$MAXIMUM_INCREM	= 2;	
	#} else {
	#		$MAXIMUM_INCREM = 1;
	#}
	
	$MAXIMUM_INCREM = ceil(abs($TCa_increm)/5.0);
	
	if ($MAXIMUM_INCREM == 0){
		$MAXIMUM_INCREM = 1;
	}

	print "Number of increments is ".$MAXIMUM_INCREM."\n"
	$INCREM = $TCa_increm/$MAXIMUM_INCREM;
	print "\nTCa increm ".$ INCREM."\n";
	$TCa = $TCa_previous
	for $i (1.. $MAXIMUM_INCREM)
	{
			fem define active;r;TCa_previous region $WALL;
			fem define initial;r;CurrentContracted region $WALL;
			fem define mapping;r;LV_CubicMapAll region $WALL;
			fem define solve;r;LV_Cubic region $WALL;

			$TCa = $TCa + $INCREM;
			MySubs::changeCalciumEspen_TCa($TCa, $i);
			$TCA_FILE = "ActivationStep/TCa_".${i};
			fem define active;r;$TCA_FILE region $WALL;
			
			print "\nTCa_goal ".${TCa_goal}."\n"
			print "TCa_previous ".${TCa_previous}."\n"
			print "TCa_Current ".${TCa}."\n"
			fem solve increment 0.0 error 1e-6 iteration 10 region $WALL;
			fem define initial;w;CurrentContracted region $WALL;
	}

	################# Output ######################
	$NAME="LV_Contracted_OptTCa";
	$FILE=$output.${NAME};
	fem define initial;w;$FILE region $WALL;

	# Export the model first
	fem export nodes;$FILE field as LVContracted_OptTCa region $WALL;
	fem export elements;$FILE field as LVContracted_OptTCa region $WALL;
	## Update the fibre field
	fem update gauss deformed_fibres region $WALL;
	fem export gauss;${FILE}."_gauss_Fibre" yg as gauss_fibre region $WALL;

	######## Save the strains and stresses to output folder ################
	fem update gauss strain extension_ratios region $WALL;
	fem export gauss;${FILE}."_gauss_ER" yg as gauss_strain region $WALL;
	fem update gauss strain region $WALL;
	fem export gauss;${FILE}."_gauss_strain" yg as gauss_strain region $WALL;
	fem update gauss strain wall region $WALL;
	fem export gauss;${FILE}."_gauss_wallstrain" yg as gauss_wallstrain region $WALL;
	fem update gauss stress total cauchy region $WALL;
	fem export gauss;${FILE}."_gauss_stress" yg as gauss_stress region $WALL;
	fem update gauss stress passive cauchy region $WALL;
	fem export gauss;${FILE}."_passive_gauss_stress" yg as gauss_stress region $WALL;
 	fem update gauss stress active cauchy region $WALL;
 	fem export gauss;${FILE}."_active_gauss_stress" yg as gauss_stress region $WALL;
	

	# Export the model again into outputAll folder. 
	$FILE2="outputAll/LVModel";
	fem export nodes;$FILE2 field as LVModel region $WALL;
	fem export eleme;$FILE2 field as LVModel region $WALL;

	# Update the cavity volume
	read commands;Cal_CavityVolumeMF; 

	
	# Write out the current pressure volume and activation level to a text file
	MySubs::writeCalciumSystole($TCa,$LV_P,$LV_CAVITY_VOLUME);
	
	##############################################################

	########### Calculate the objective function ###################
	read commands;CalObjective;

	fem quit


	
