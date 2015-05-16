use POSIX

##### Auxiliary parameters
$ITERS=20;   # Newton-Raphson iteration number. 
$ERROR_TOLERANCE=1e-6;

# Read in the previous solved solution
fem define acti;r;TCa_previous region $WALL;
fem define initial;r;CurrentContracted region $WALL;
fem evaluate residual;CurrentContracted wrt geom_params region $WALL;
system('tail CurrentContracted.opresi');

# fem solve increm 0.0 iter 10 error 1e-6 region $WALL;
# Use pressure increment and displacement boundary conditions in ipinit file to run simulation. 
fem define initial;r;LV_CubicPreEpiBase region $WALL;
fem define mapping;r;LV_CubicMapAll region $WALL;
fem define solve;r;LV_Cubic region $WALL;

set echo
####### Extract total pressure load. ########################################
$Start_DS=`sed -e s%D%E% LV_CubicPreEpiBase.ipinit | awk -v line=16 'NR==line{printf("%d",\$1)}'`;
if ($Start_DS==1)
{
	$LV_P_goal=`sed -e s%D%E% LV_CubicPreEpiBase.ipinit | awk -v line=360 'NR==line{printf("%.5f",\$5)}'`;
	$TCa_goal=`sed -e s%D%E% TCa_current.ipacti | awk -v line=16 'NR==line{printf("%.5f",\$7)}'`;
	$TCa_previous = `sed -e s%D%E% TCa_previous.ipacti | awk -v line=16 'NR==line{printf("%.5f",\$7)}'`;
} else {
	$LV_P_goal=`sed -e s%D%E% LV_CubicPreEpiBase.ipinit | awk -v line=531 'NR==line{printf("%.5f",\$5)}'`;
	$TCa_goal=`sed -e s%D%E% TCa_current.ipacti | awk -v line=16 'NR==line{printf("%.5f",\$7)}'`;
	$TCa_previous = `sed -e s%D%E% TCa_previous.ipacti | awk -v line=16 'NR==line{printf("%.5f",\$7)}'`;
}

$TCa_step = $TCa_goal - $TCa_previous;

print "\033[0;30;45m Solving from TCa = $TCa_previous to reach TCa = $TCa_goal, with increment = $TCa_step    \033[0m\n";
print "\033[0;30;45m Solving to apply pressure = $LV_P_goal    \033[0m\n";

####### Adapt number of increments to size of total pressure load step ##################
$MAXIMUM_INCREM = ceil(abs($TCa_step)/2.0);
	
if ($MAXIMUM_INCREM == 0){
	$MAXIMUM_INCREM = 1;
}
print $MAXIMUM_INCREM."\n"
$TCa = $TCa_previous
###### Start for loop to increment pressure by reasonable intermediate steps ############
for $i ( 1..$MAXIMUM_INCREM) 
{
    $TCa_increm = ($TCa_step)/$MAXIMUM_INCREM;
    # Apply TCa increment
    $TCa= $TCa + $TCa_increm;
    MySubs::changeCalcium_TCa($TCa, $i);
	
    $TCA_FILE = "output_debug/TCa_".${i};
    fem define active;r;$TCA_FILE region $WALL;
    print "\033[0;30;42m    Goal: TCa = $TCa_goal kPa    \033[0m\n";
    print "\033[0;30;43m    Apply: TCa increm $TCa_increm kPa    \033[0m\n";
    print "\033[0;30;47m    Load step: $i out of $MAXIMUM_INCREM    \033[0m\n";   
    fem solve increment 0.0 iter $ITERS error $ERROR_TOLERANCE region $WALL;
    print "\033[0;30;46m    Current: TCa = $TCa kPa    \033[0m\n";
    read commands;Extract_Pressure;

    $P_INCREM = 1/$MAXIMUM_INCREM
    # Solve with defined increment to apply desired boundary conditions
    read commands;Extract_Pressure;
    $LVP_increm = $LV_P_goal*$P_INCREM;
    print "\033[0;30;42m    Goal: apply pressure = $LV_P_goal kPa    \033[0m\n";
    print "\033[0;30;43m    Apply: pressure increm $LVP_increm kPa    \033[0m\n";
    print "\033[0;30;47m    Load step: $i out of $MAXIMUM_INCREM    \033[0m\n";
    fem solve increment $P_INCREM iter $ITERS  error $ERROR_TOLERANCE;
    read commands;Extract_Pressure;

    ########################## Output ###############################
    $NAME="LVSystole_".${outputDebug_index};
    $FILE=$outputDebug.${NAME};

    # Save the current status as the initial condition for the next phase
    fem define initial;w;$FILE region $WALL;

    # Export the model into output folder.
    fem export nodes;$FILE field as LVSystole region $WALL;
    fem export elements;$FILE field as LVSystole region $WALL;

    ######## Save the strains and stresses to output folder ################
    fem update gauss strain extension_ratios region $WALL;
    fem export gauss;${FILE}."_gauss_ER" yg as gauss_strain region $WALL;
    fem update gauss strain region $WALL;
    fem export gauss;${FILE}."_gauss_strain" yg as gauss_strain region $WALL;
    fem update gauss strain wall region $WALL;
    fem export gauss;${FILE}."_gauss_wall_strain" yg as gauss_wall_strain region $WALL;
    fem update gauss stress total cauchy region $WALL;
    fem export gauss;${FILE}."_gauss_stress" yg as gauss_stress region $WALL;
    fem update gauss stress passive cauchy region $WALL;
    fem export gauss;${FILE}."_passive_gauss_stress" yg as gauss_stress region $WALL;
    fem update gauss stress active cauchy region $WALL;
    fem export gauss;${FILE}."_active_gauss_stress" yg as gauss_stress region $WALL;

    $outputDebug_index = $outputDebug_index + 1;

    ############### Redefine initial condition and insolv file for the wall
    fem define init;r;$FILE region $WALL;
    fem define mapping;r;LV_CubicMapAll region $WALL;
    fem define solve;r;LV_Cubic region $WALL;

    $str = "echo ".${i}." > loaditeration.txt"
    system($str)		
}

$str = "echo ".${CONVERGED}." > convergence.txt"
#$str = "echo 1 > convergence.txt"
print $str
system($str)

MySubs::changeCalcium_TCa($TCa, "current");
#fem define acti;w;${output}."TCa" region $WALL;
fem export nodes;${output}."LVContraction" field as LVInflation region $WALL;
fem export elements;${output}."LVContraction" field as LVInflation region $WALL;
fem define initial;w;${output}."LVContraction" region $WALL

#Update extension ratio and strain
fem update gauss strain extension_ratios region $WALL;
fem export gauss;${output}."LVContraction"."_gauss_ER" yg as gauss_strain region $WALL;
fem update gauss strain region $WALL;
fem export gauss;${output}."LVContraction"."_gauss_strain" yg as gauss_strain region $WALL;
fem update gauss stress total cauchy region $WALL;
fem export gauss;${output}."LVContraction"."_gauss_stress" yg as gauss_stress region $WALL;
fem update gauss stress passive cauchy region $WALL;
fem export gauss;${output}."LVContraction"."_passive_gauss_stress" yg as gauss_stress region $WALL;
fem update gauss stress active cauchy region $WALL;
fem export gauss;${output}."LVContraction"."_active_gauss_stress" yg as gauss_stress region $WALL;
fem update gauss deformed_fibres region $WALL;
fem export gauss;${output}."LVContraction"."_gauss_Fibre" yg as gauss_fibre region $WALL;

fem list stress;${output}."gauss_stress" full region $WALL;
fem list strain;${output}."gauss_strain" full region $WALL;
fem list stress;${output}."passive_gauss_stress" passive full region $WALL;
fem list stress;${output}."active_gauss_stress" active full region $WALL;

	
################## Calculate the cavity volume ########################
read commands;Extract_Pressure
read commands;Cal_CavityVolume
