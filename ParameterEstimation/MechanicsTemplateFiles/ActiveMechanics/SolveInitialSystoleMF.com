### For python monitoring purposes, set the finished boolean as false in the beginning. 
$str = "echo 0 > loaditeration.txt"
system($str)

# Default converged parameter should be zero. 
$str = "echo 0 > convergence.txt"
system($str)

############ Set up output directory ##########################
read commands;SetOutput;

########### Read in the reference wall model #####################
read commands;ReadRefWallModel;

########### Read in the reference cavity model ###################
read commands;ReadRefCavityModel;

read commands;ReadInflated;
read commands;SolveIVCEjectionMF;

fem list elem;output_cavity_volume/LVCavityCurrent deformed total region $LV_CAVITY

############ Project registered surface data #######################
read commands;CalObjective;

### Once file has finished executing all commands, set finished boolean to true. 
$str = "echo 1 > finished.txt"
system($str)

fem quit;
