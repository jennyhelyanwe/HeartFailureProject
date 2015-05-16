########Must do this command to ensure that the projection is to the deformES surface ########
fem update geometry from solution region $WALL;
	
$FILE = ${outputDebug}."ES"
fem exp node;$FILE as ES region $WALL;
fem exp elem;$FILE as ES region $WALL;

fem def nodes;w;$FILE region $WALL;
fem def elem;w;$FILE region $WALL;

read commands;CreateSurfacePoints

system('python RegisterData.py')

# Rigidly transform the surface data
CORE::open( MYFILE,'output_debug/TransformationMatrix.TRN' ) or die( "ERROR: Unable to open output file\n" );
while (<MYFILE>) {
 	chomp;
 	$transformation= "$_\n";
 }
close (MYFILE); 
@temp = split(',',$transformation);
$translate_x=@temp[0]
$translate_y=@temp[1]
$translate_z=@temp[2]
$rotation_x=@temp[3]
$rotation_y=@temp[4]
$rotation_z=@temp[5]
$pi = 3.1415926
$rotation_x_deg = $rotation_x*180/$pi
$rotation_y_deg = $rotation_y*180/$pi
$rotation_z_deg = $rotation_z*180/$pi

$translate=$translate_x.",".$translate_y.",".$translate_z
$rotation = $rotation_x_deg.",".$rotation_y_deg.",".$rotation_z_deg
print $translate
print $rotation

##### Rotate and translate the epi data #######################
fem def data;r;Surface_Points_Epi_ES region $WALL;
fem exp data;${outputDebug}."ES_Epi_nonreg" as nonreg region $WALL;

fem change data rotate by $rotation_z_deg,0,0 about 0,0,0 axis 0,0,1;
fem change data rotate by $rotation_y_deg,0,0 about 0,0,0 axis 0,1,0;
fem change data rotate by $rotation_x_deg,0,0 about 0,0,0 axis 1,0,0;

fem change data translate by $translate_x,$translate_y,$translate_z
fem exp data;${outputDebug}."ES_Epi_reg" as trans region $WALL;
fem def data;w;${outputDebug}."Surface_Points_Epi_reg" region $WALL;

###### Rotate and translate the endo data ###################
fem def data;r;Surface_Points_Endo_ES region $WALL;
fem exp data;${outputDebug}."ES_Endo_nonreg" as nonreg region $WALL;

fem change data rotate by $rotation_z_deg,0,0 about 0,0,0 axis 0,0,1;
fem change data rotate by $rotation_y_deg,0,0 about 0,0,0 axis 0,1,0;
fem change data rotate by $rotation_x_deg,0,0 about 0,0,0 axis 1,0,0;

fem change data translate by $translate_x,$translate_y,$translate_z
fem exp data;${outputDebug}."ES_Endo_reg" as trans region $WALL;
fem def data;w;${outputDebug}."Surface_Points_Endo_reg" region $WALL;

######### Evaluating Epi Error ###################################
fem define data;r;${outputDebug}."Surface_Points_Epi_reg" region $WALL;
fem group faces 5,9,13,16,21,25,29,32,37,41,45,48,53,57,61,64 as EPI region $WALL;
fem def xi;c closest_face faces EPI search_start 20 region $WALL;

fem list data;${error}."EpiProjectionToES" error
fem exp data;${error}."EpiProjectionToES" as EpiProjectionToES error region $WALL;

######### Evaluating Endo error #############################
fem def data;r;${outputDebug}."Surface_Points_Endo_reg" region $WALL;
fem group faces 4,8,12,15,20,24,28,31,36,40,44,47,52,56,60,63 as ENDO region $WALL;
fem def xi;c closest_face faces ENDO region $WALL;

fem list data;${error}."EndoProjectionToES" error
fem export data;${error}."EndoProjectionToES" as EndoProjectionToES error region $WALL;

