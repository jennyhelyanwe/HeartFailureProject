use POSIX;
set echo on;

############ Set up output directory ##########################
read commands;SetOutput;

########### Read in the reference wall model #####################
read commands;ReadRefWallModel;
	
########### Read in the reference cavity model ###################
read commands;ReadRefCavityModel;
	
#################### Decide on how many steps to take ################
$C1_current = `sed -e s%D%E% LV_CubicGuc.ipmate | awk -v line=41 'NR==line{printf("%.5f",\$5)}'` 
$C1_previous = `sed -e s%D%E% LV_CubicGuc_previous.ipmate | awk -v line=41 'NR==line{printf("%.5f",\$5)}'` 
        
print "\nPrevious C1 is ".${C1_previous}."\n"
print "C1 goal is ".${C1_current}."\n"
$C1_incr = $C1_current - $C1_previous;

print "\n \033[0;30;42m ================================================================================== \033[0m\n";
print "\033[0;30;42m       Solving iteratively from c1 = $C1_previous to reach c1 = $C1_current    	\033[0m\n";
print "\033[0;30;42m ================================================================================== \033[0m\n";

$no_steps = ceil(abs($C1_incr)/1.0)
	
if ($no_steps ==0) {
	$no_steps = 1;
}
print "\n*********** Number of steps is ".${no_steps}."\n"
$C1_subincr = $C1_incr/$no_steps;
print "\n*********** The C1_smallincr is \n"
print $C1_subincr
	
$C1 = $C1_previous
print "\nCurrent starting C1 is ".${C1}."\n"
for ($i = 1; $i <= $no_steps; $i++){
	fem define initial;r;CurrentInflated region $WALL;
	fem define mapping;r;LV_CubicMapAll region $WALL;  
	fem define solve;r;LV_Cubic region $WALL;

	$C1 = $C1 + $C1_subincr;
	print "\nCurrent C1 is equal to: ".${C1}."\n"
	# Rewrite ipmate file
	$cmd = "python createIPMATE.py ".${i}." '".${C1}."'"
	system($cmd)
		
	# Re-read in ipmate file
	$FILE = ${outputDebug}."LV_CubicGuc_".${i}
	fem define mate;r;$FILE region $WALL;
	fem solve increment 0.0 error 1e-6 iteration 10 region $WALL;
	fem define initial;w;CurrentInflated region $WALL;
}

################# Output ######################
$NAME="LVInflation";
$FILE=$output.${NAME};
# Export the model first
fem export nodes;$FILE field as LVInflation_Update region $WALL;
fem export elements;$FILE field as LVInflation_Update region $WALL;
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

fem list stress;${output}."gauss_stress" full region $WALL;
fem list strain;${output}."gauss_strain" full region $WALL;
fem list stress;${output}."passive_gauss_stress" passive full region $WALL;

# Update the cavity volume
read commands;Extract_Pressure;
read commands;Cal_CavityVolume;

fem list elem;output_cavity_volume/LVCavityUpdate deformed total region $LV_CAVITY

########### Calculate the objective function ###################
read commands;CalObjective;

fem quit;


	
