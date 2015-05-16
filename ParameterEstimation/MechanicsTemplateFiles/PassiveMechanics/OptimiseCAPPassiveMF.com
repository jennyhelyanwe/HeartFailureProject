########################################################################################
################################# Optimization #########################################
########################################################################################
use POSIX
# This program is designed to optimize the stiffness parameter C1 of the transversely-isotropic
# law using a set of 3D material points obtained from tag tracking. The target data points will
# be material points tracked at frame 6 (hence end-diastolic state), and the landmark points 
# are points from frame 48. 

########################### Set up initial files ######################################



	set echo on;
	

	############ Set up output directory ##########################
	read commands;SetOutput;

	########### Read in the reference wall model #####################
	read commands;ReadRefWallModel;
	
	########### Read in the reference cavity model ###################
	read commands;ReadRefCavityModel;
	
	
	#################### Read in pre-solved ED model ################

	
	
	system('cp LV_CubicGuc_temp.ipmate LV_CubicGuc.ipmate');

	$C1_current = `sed -e s%D%E% LV_CubicGuc.ipmate | awk -v line=41 'NR==line{printf("%.5f",\$5)}'` 
	$C1_previous = `sed -e s%D%E% LV_CubicGuc_previous.ipmate | awk -v line=41 'NR==line{printf("%.5f",\$5)}'` 
        
	print "\nPrevious C1 is ".${C1_previous}."\n"
	print "C1 goal is ".${C1_current}."\n"
	$C1_incr = $C1_current - $C1_previous;
	#if (abs($C1_incr)>0.5 && abs($C1_incr) < 1){
    #    	$no_steps = 2;
	#} elsif (abs($C1_incr)>=1&& abs($C1_incr)< 2){
	#	$no_steps = 3;
    #} elsif (abs($C1_incr)>=2){
    #    $no_steps = 5;
	#} else {
	#	$no_steps = 1;
	#}

	$no_steps = ceil(abs($C1_incr)/0.3)
	
	if ($no_steps ==0) {
		$no_steps = 1;
	}
	print "\n*********** Number of steps is ".${no_steps}."\n"
	$C1_smallincr = $C1_incr/$no_steps;
	print "\n*********** The C1_smallincr is \n"
	print $C1_smallincr
	
	$C1 = $C1_previous
	print "\nCurrent starting C1 is ".${C1}."\n"

	for ($i = 1; $i <= $no_steps; $i++){
		fem define initial;r;CurrentInflated region $WALL;
		fem define mapping;r;LV_CubicMapAll region $WALL;  
		fem define solve;r;LV_Cubic region $WALL;

		$C1 = $C1 + $C1_smallincr;
		print "\nCurrent C1 is equal to: ".${C1}."\n"
		# Rewrite ipmate file
		$cmd = "python createIPMATE.py ".${i}." '".${C1}."'"
		system($cmd)
		
		# Re-read in ipmate file
		$FILE = "LV_CubicGuc_".${i}
		fem define mate;r;$FILE region $WALL;
		fem solve increment 0.0 error 1e-6 iteration 10 region $WALL;

		fem define initial;w;CurentInflated region $WALL;

	}
	
	fem define initial;w;CurrentInflatedUpdated region $WALL;

	################# Output ######################
	$NAME="LV_Inflation_OptC1";
   	$FILE=$output.${NAME};
    	fem define initial;w;$FILE region $WALL;

	# Export the model first
	fem export nodes;$FILE field as LVInflation_OptC1 region $WALL;
    	fem export elements;$FILE field as LVInflation_OptC1 region $WALL;
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

	$FILE2="outputAll/LVModel";
	fem export nodes;$FILE2 field as LVModel region $WALL;
	fem export eleme;$FILE2 field as LVModel region $WALL;

	# Update the cavity volume
        read commands;Extract_Pressure;
	read commands;Cal_CavityVolume;

	fem list elem;outputCavityVolume/LVCavity_Current deformed total region $LV_CAVITY
	
	##############################################################

	########### Calculate the objective function ###################
	read commands;CalObjective;

	fem quit


	
