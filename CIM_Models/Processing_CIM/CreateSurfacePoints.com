
##################################### ED #####################################
#### Read in the set-up files
fem define para;r;small;			# Define the parameter file
fem define coor;r;mapping;			# Define the coordinate file which involves mapping
fem define base;r;LVBasisV2

#### Read in the geometric model
fem define node;r;fitted_endo_Phase48_CIM_CorrectApex_UpS3
fem define elem;r;fitted_endo_Phase48_CIM_CorrectApex_UpS3


#### Read in the local FE coordinates for endocardial surface #######
fem define data;r;Temp			## This is just a temporary data file
fem define xi;r;Surface_Points_Endo
fem define data;c from_xi
fem define data;w;Surface_Points_Endo_ED
fem list data statistics
fem export data;Surface_Points_Endo_ED as Surface_Points_Endo_ED


#### Read in the local FE coordinates for epicardial surface #######
fem define data;r;Temp			## This is just a temporary data file
fem define xi;r;Surface_Points_Epi
fem define data;c from_xi
fem define data;w;Surface_Points_Epi_ED
fem list data statistics
fem export data;Surface_Points_Epi_ED as Surface_Points_Epi_ED


fem reallocate

##################################### ES #####################################
#### Read in the set-up files
fem define para;r;small;			# Define the parameter file
fem define coor;r;mapping;			# Define the coordinate file which involves mapping
fem define base;r;LVBasisV2

#### Read in the geometric model
fem define node;r;fitted_endo_Phase48_CIM_CorrectApex_UpS3
fem define elem;r;fitted_endo_Phase48_CIM_CorrectApex_UpS3


#### Read in the local FE coordinates for endocardial surface #######
fem define data;r;Temp			## This is just a temporary data file
fem define xi;r;Surface_Points_Endo
fem define data;c from_xi
fem define data;w;Surface_Points_Endo_ES
fem list data statistics
fem export data;Surface_Points_Endo_ES as Surface_Points_Endo_ES


#### Read in the local FE coordinates for epicardial surface #######
fem define data;r;Temp			## This is just a temporary data file
fem define xi;r;Surface_Points_Epi
fem define data;c from_xi
fem define data;w;Surface_Points_Epi_ES
fem list data statistics
fem export data;Surface_Points_Epi_ES as Surface_Points_Epi_ES


fem quit


