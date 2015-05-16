########## Read in the reference cavity model ##############

########## Define the Cavity Mesh ###########
fem define region;r;LV_Cubic  
fem define coordinates;r;LV_Cavity region $LV_CAVITY			# Define coordinate system
                      
fem define nodes;r;DSCavity region $LV_CAVITY		# Define the Cavity mesh
fem define elements;r;DSCavity region $LV_CAVITY		# Define the Cavity mesh

# Read in the scale factors of the LV cavity
fem define base;r;LV_CubicReadSF
fem define line;r;DSWall
fem define base;r;LV_Cubic
	
fem export node;DSCavity  as ReferenceCavityModel region $LV_CAVITY
fem export element;DSCavity as ReferenceCavityModel region  $LV_CAVITY

fem list element total region $LV_CAVITY;
fem list element;DSCavity total region $LV_CAVITY;
	
# Define the fibre field
fem define fibre;r;DTMRI_CIMModel region $WALL;		
fem define elem;r;DTMRI_CIMModel fibre region $WALL;

fem export node;UPFInitialModelWithFibre as UPFInitialModel region $WALL;
fem export elem;UPFInitialModelWithFibre as UPFInitialModel region $WALL;

# Define the equation file
fem define equation;r;LV_Cubic region $WALL;
# Define the material file
fem define material;r;LV_CubicGuc region $WALL;
# Definie active file
fem define active;r;LV_CubicCaActivIni region $WALL;
	
fem define equation;r;cavity region $LV_CAVITY
fem define initial;r;cavity region $LV_CAVITY
