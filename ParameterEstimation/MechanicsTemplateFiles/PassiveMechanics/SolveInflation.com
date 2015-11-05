use POSIX;
################################# Solve Passive Inflation ################

set echo;
	#==================================================================================================
print "\033[0;30;42m ============================================================ \033[0m\n"; 
print "\033[0;30;42m                  Initial Conditions                          \033[0m\n"; 
print "\033[0;30;42m ============================================================ \033[0m\n"; 
	#==================================================================================================

fem define initial;r;LV_CubicPreEpiBase region $WALL;

	#==================================================================================================
print "\033[0;30;42m ============================================================ \033[0m\n"; 
print "\033[0;30;42m                    Solution Method                           \033[0m\n"; 
print "\033[0;30;42m ============================================================ \033[0m\n"; 
	#==================================================================================================


fem define mapping;r;LV_CubicMapAll region $WALL;
fem define solve;r;LV_Cubic region $WALL;  

########### Check the system is converged before solving anything
set output;InflatioOutput on

fem list variables;pressure/Pressure_Initial region $WALL;

########## List the LV cavity volume at the reference state
fem list elements;output_cavity_volume/LVCavityInit total region $LV_CAVITY; # Absolutely is correct
$LV_CAVITY_VOLUME=`awk '\$1=="Total"{printf("%f",\$4);exit}\' output_cavity_volume/LVCavityInit.opelem`; 
print "\n\033[0;30;42mDS LV cavity volume = ${LV_CAVITY_VOLUME}\033[0m\n";

#==========================================================================================================
print "\033[0;30;42m ============================================================ \033[0m\n"; 
print "\033[0;30;42m                       Solve Inflation                        \033[0m\n"; 
print "\033[0;30;42m ============================================================ \033[0m\n"; 
#===================================================================================
print "\033[0;30;46m Increase cavity pressures incrementally to simulate diastole \033[0m\n"; 

   	
$output_index=0;
$NAME="LVInflation_".${output_index};
$FILE=$outputDebug.${NAME};

# Export the model first
fem export nodes;$FILE as LVInflation field region $WALL;
fem export elements;$FILE as LVInflation field region $WALL;
$output_index=$output_index+1;

$output_index=1;
$ITERS =20;
$TOL = 1e-6;

####### Check whether solving starts from DS ##########
$Start_DS=`sed -e s%D%E% LV_CubicPreEpiBase.ipinit | awk -v line=16 'NR==line{printf("%d",\$1)}'`;
if ($Start_DS==1)
{
	$LV_P=`sed -e s%D%E% LV_CubicPreEpiBase.ipinit | awk -v line=360 'NR==line{printf("%.5f",\$5)}'`;
} else {
	$LV_P=`sed -e s%D%E% LV_CubicPreEpiBase.ipinit | awk -v line=531 'NR==line{printf("%.5f",\$5)}'`;
}
print "\033[0;30;45m Applying additional pressure = $LV_P   \033[0m\n";

$MAXIMUM_INCREM = ceil(abs($LV_P)/0.05)
if ($MAXIMUM_INCREM==0.0)
{
	$MAXIMUM_INCREM= 1
}

if ($MAXIMUM_INCREM==1)
{
    $MAXIMUM_INCREM=2
}

if ($LV_P < 0)
{
	$MAXIMUM_INCREM = $MAXIMUM_INCREM + 1
}

print "\033[0;30;42mThis is going to take ".$MAXIMUM_INCREM." load step(s).\033[0m\n"
$INCREM = 1/$MAXIMUM_INCREM;
$p_increm = $INCREM*$LV_P
for $i ( 1..$MAXIMUM_INCREM ) 
{       
    print "\n\033[0;30;42mSolve load step number ".$i."\033[0m\n"
	
	print "\033[0;30;43m    Apply: pressure increment $p_increm kPa    \033[0m\n";
    print "\033[0;30;47m    Load step: $i out of $MAXIMUM_INCREM    \033[0m\n";
	fem solve increment $INCREM iter $ITERS error $TOL region $WALL;
       	 	
	fem list variables;Pressure_$i region $WALL;

	read commands;Extract_Pressure;
	$FILE = ${NAME1}.${LV_P};

	$NAME="LVInflation_".${output_index};
	$FILE=$outputDebug.${NAME};
	fem define initial;w;$FILE region $WALL;
	
	# Export the model first
	fem export nodes;$FILE field as LVInflation region $WALL;
	fem export elements;$FILE field as LVInflation region $WALL;

	$output_index=$output_index+1;

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

	############### Redefine initial condition and insolv file for the wall
	fem define initial;r;$FILE region $WALL;
	fem define mapping;r;LV_CubicMapAll region $WALL;  
	fem define solve;r;LV_Cubic region $WALL;       	
}

fem export nodes;${output}."LVInflation" field as LVInflation region $WALL;
fem export elements;${output}."LVInflation" field as LVInflation region $WALL;

fem define initial;w;${output}."LVInflation" region $WALL

# Update extension ratio and strain
fem update gauss strain extension_ratios region $WALL;
fem export gauss;${output}."LVInflation"."_gauss_ER" yg as gauss_strain region $WALL;
fem update gauss strain region $WALL;
fem export gauss;${output}."LVInflation"."_gauss_strain" yg as gauss_strain region $WALL;
fem update gauss stress total cauchy region $WALL;
fem export gauss;${output}."LVInflation"."_gauss_stress" yg as gauss_stress region $WALL;
fem update gauss stress passive cauchy region $WALL;
fem export gauss;${output}."LVInflation"."_passive_gauss_stress" yg as gauss_stress region $WALL;
fem update gauss stress active cauchy region $WALL;
fem export gauss;${output}."LVInflation"."_active_gauss_stress" yg as gauss_stress region $WALL;
fem update gauss deformed_fibres region $WALL;
fem export gauss;${output}."LVInflation"."_gauss_Fibre" yg as gauss_fibre region $WALL;

fem list stress;${output}."gauss_stress" full region $WALL;
fem list strain;${output}."gauss_strain" full region $WALL;
fem list stress;${output}."passive_gauss_stress" passive full region $WALL;

################## Calculate the cavity volume ########################
read commands;Extract_Pressure
read commands;Cal_CavityVolume
MySubs::writePVdiastole($LV_P,$LV_CAVITY_VOLUME);

