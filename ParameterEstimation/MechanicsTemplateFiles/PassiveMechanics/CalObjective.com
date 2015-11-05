############## Must do this command to ensure that the projection is to the deformed surface ##########	
fem update geometry from solution region $WALL;
	
		
############ Define landmark and target data points ###############
##################################### Epi ##############################
# Define the undata (landmark points)
fem define data;r;Surface_Points_Epi_ED region $WALL;	
fem exp data;${outputDebug}."ED_Epi_nonreg" as nonreg region $WALL;

fem export data;${error}."ED_Epi" as ED_Epi region $WALL;	
#fem list data statistics
	
## Project ED Epi data to epi surface first
fem group faces 5,9,13,16,21,25,29,32,37,41,45,48,53,57,61,64 as EPI region $WALL;
fem def xi;c closest_face faces EPI search_start 20 region $WALL;
	
## List the projection error
fem list data;${error}."EpiProjectionToED" error
fem export data;${error}."EpiProjectionToED" as EpiProjectionToED error region $WALL;
	
	
##################################### Endo ##############################
# Define the undata (landmark points)
fem define data;r;Surface_Points_Endo_ED region $WALL;	
fem exp data;${outputDebug}."ED_Endo_nonreg" as nonreg region $WALL;

fem export data;${error}."ED_Endo" as ED_Endo region $WALL;	
#fem list data statistics
	
## Project CAP ED Endo data to endo surface first
fem group faces 4,8,12,15,20,24,28,31,36,40,44,47,52,56,60,63 as ENDO region $WALL;
	
fem def xi;c closest_face faces ENDO search_start 20  region $WALL;
	
## List the projection error
fem list data;${error}."EndoProjectionToED" error
fem export data;${error}."EndoProjectionToED" as EndoProjectionToED error region $WALL;

