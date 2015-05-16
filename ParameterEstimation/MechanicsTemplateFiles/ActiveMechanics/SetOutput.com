	
use MySubs;
#################################### Output Directory ##############################################
# Set up the output directory for storing output files
	
$output = "output/";
if( ! -d ${output})
{
  mkdir ${output};
}

$outputDebug = "output_debug/";
if( ! -d ${outputDebug})
{
  mkdir ${outputDebug};
}

$outputCavityVolume="output_cavity_volume/";
if( ! -d ${outputCavityVolume})
{
    mkdir ${outputCavityVolume};
}	

$Pressure="pressure/";
if( ! -d ${Pressure})
{
    mkdir ${Pressure};
}

$error = "output_errors/";
if( ! -d ${error})
{
    mkdir ${error};
}
	
# Initialize some logical variables

$TRUE = 1;  
$FALSE = 0;
	
#$SOLVE_INFLATION = $TRUE;
$SOLVE_INFLATION = $FALSE; 					# Enable if simulation is not required

#$SOLVE_IVC = $TRUE;
$SOLVE_IVC = $FALSE;
	
$SOLVE_EJECTION = $TRUE;						# Enable if simulation is not required
#$SOLVE_EJECTION = $FALSE;
	
#$SOLVE_IVR = $TRUE;						# Enable if simulation is not required
#$SOLVE_IVR = $FALSE;
	
#$count=1;
