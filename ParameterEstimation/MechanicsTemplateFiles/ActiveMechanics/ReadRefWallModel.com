	
########## Read in the reference model ############
# Define the parameters and coordinate and regions
fem define para;r;LV_Cubic
fem define coor;r;LV_Cubic
fem define region;r;LV_Cubic

# Define the regions
$WALL=1;
$LV_CAVITY=2;
	
# Define the basis functions
fem define base;r;LV_Cubic
	
# Define the geometry of the model
fem def node;r;DSWall region $WALL
fem def elem;r;DSWall region $WALL
	
# Define the basis functions which allow the scale factors to be read in 
fem define base;r;LV_CubicReadSF
fem define line;w;DSWall
fem define base;r;LV_Cubic
	
fem list element total region $WALL;
# Output the reference model
fem export node;DSWall as ReferenceWallModel region $WALL;	
fem export element;DSWall as ReferenceWallModel region $WALL;	
	
########## Set up Nodal Groups ############
# Some nodes have been grouped based on different purposes such as fixed position, nodal group where particular boundary conditions will be applied, nodal group which exhibit different material properties. 

fem group nodes 14..17 as ENDOBASENODES;
fem group nodes 31..34 as EPIBASENODES;
fem group nodes 14..17,31..34 as BASENODESALL;
	
########## Set up element groups ###########
# Elements were also grouped based on different regions 
fem group elements 1..16 as WALL region $WALL;

	
