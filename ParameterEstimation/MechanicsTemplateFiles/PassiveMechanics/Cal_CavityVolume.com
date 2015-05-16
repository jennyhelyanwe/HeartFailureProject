################################################# Compute LV Volume ########################################################

fem define coupling;r;coupled

fem define solve;r;coupled coupled region $WALL,$LV_CAVITY; # Define solve for coupled regions	

fem update solution coupled source_region $WALL; # Copy deformed endocardial surfaces (interface nodes) from wall region to cavity region
		
# Update locations of cavity nodes based on average positions of endocardial 'rings' (in 1 represents the RC coordinate along which node 118 is averaged among nodes 14..17)

fem update solution cavity_reference average 118 in 1 node 14..17 region $LV_CAVITY;
			
####################################### Calculates LV Cavity Volume #########################################
fem list elements;output_cavity_volume/LVCavity_${LV_P} deformed total region $LV_CAVITY;
$LV_CAVITY_VOLUME=`awk '\$1=="Total"{printf("%d",\$4);exit}\' output_cavity_volume/LVCavity_${LV_P}.opelem`;
print "\n\033[0;30;42mLV pressure = ${LV_P} LV_volume=${LV_CAVITY_VOLUME} \033[0m\n";

fem exp node;output_cavity_volume/LVCavity_${LV_P} as LVInflation field region $LV_CAVITY;
fem exp elem;output_cavity_volume/LVCavity_${LV_P} as LVInflation field region $LV_CAVITY;

########## Redefine files required for the wall mechanics simulation ############################
fem define coupling;r;coupled_reset


