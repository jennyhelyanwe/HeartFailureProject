set echo on;
system('rm *.*~');

###################################### Define auxiliary parameters ################################

###################################################################################################

####################################### FEM starts ################################################
fem def para;r;small;			# Define the parameter file
fem def coor;r;mapping;			# Define the coordinate file which involves mapping

##################### Read in the model at DS ################################
fem define node;r;DSWall
fem def base;r;LVBasisV2;			# Define the basis functions
fem define elem;r;DSWall

## Output the DS model
fem export node;DSWall as DSWall;
fem export elem;DSWall as DSWall;

##################### Read in the cavity model #########################
fem define node;r;DSCavity
fem define elem;r;DSCavity

fem export node;DSCavity as DSCavity
fem export elem;DSCavity as DSCavity

fem quit
