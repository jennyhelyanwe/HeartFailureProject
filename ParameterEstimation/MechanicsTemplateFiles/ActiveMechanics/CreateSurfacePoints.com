#### This com file creates surface data points for the generic geometric model based on the number of element points specified by the ipxi file. 

#### Read in the set-up files
#fem def para;r;LV_Cubic;
fem define para;r;small;		# Define the parameter file
fem define coor;r;mapping;		# Define the coordinate file which involves mapping
fem define base;r;LVBasisV2

#### Read in the geometric model
$FILE = ${outputDebug}."ES"
fem define node;r;$FILE
fem define elem;r;$FILE

#### Read in the local FE coordinates for endocardial surface #######
fem define data;r;Temp			# This is just a temporary data file
fem define xi;r;Surface_Points_Endo
fem define data;c from_xi
fem define data;w;${outputDebug}."ES_model_surface_endo"
fem list data statistics
fem export data;${outputDebug}."ES_model_surface_endo" as ES_model_surface_endo

#### Read in the local FE coordinates for epicardial surface #######
fem define data;r;Temp			# This is just a temporary data file
fem define xi;r;Surface_Points_Epi
fem define data;c from_xi
fem define data;w;${outputDebug}."ES_model_surface_epi"
fem list data statistics
fem export data;${outputDebug}."ES_model_surface_epi" as ES_model_surface_epi

#fem reallocate

