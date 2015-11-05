############## Must do this command to ensure that the projection is to the deformES surface ##########	
fem update geometry from solution region $WALL;
	
		
############ Define landmark and target data points ###############
##################################### Epi ##############################
# Define the undata (landmark points)
fem define data;r;Surface_Points_Epi_ES region $WALL;	
fem exp data;${outputDebug}."ES_Epi_nonreg" as nonreg region $WALL;

fem export data;${error}."ES_Epi" as ES_Epi region $WALL;	
#fem list data statistics
	
## Project ES Epi data to epi surface first
fem group faces 5,9,13,16,21,25,29,32,37,41,45,48,53,57,61,64 as EPI region $WALL;
fem def xi;c closest_face faces EPI search_start 20 region $WALL;
	
## List the projection error
fem list data;${error}."EpiProjectionToES" error
fem export data;${error}."EpiProjectionToES" as EpiProjectionToES error region $WALL;
	
	
##################################### Endo ##############################
# Define the undata (landmark points)
fem define data;r;Surface_Points_Endo_ES region $WALL;	
fem exp data;${outputDebug}."ES_Endo_nonreg" as nonreg region $WALL;
	
fem export data;${error}."ES_Endo" as ES_Endo region $WALL;	
#fem list data statistics
	
## Project CAP ES Endo data to endo surface first
fem group faces 4,8,12,15,20,24,28,31,36,40,44,47,52,56,60,63 as ENDO region $WALL;
	
fem def xi;c closest_face faces ENDO search_start 20  region $WALL;
	
## List the projection error
fem list data;${error}."EndoProjectionToES" error
fem export data;${error}."EndoProjectionToES" as EndoProjectionToES error region $WALL;

