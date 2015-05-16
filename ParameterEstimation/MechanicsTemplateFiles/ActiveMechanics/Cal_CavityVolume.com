################################################# Compute LV Volume ########################################################
fem define equation;r;cavity region $LV_CAVITY
fem define initial;r;cavity region $LV_CAVITY

fem define coupling;r;coupled

fem define solve;r;coupled coupled region $WALL,$LV_CAVITY; # Define solve for coupled regions
fem update solution coupled source_region $WALL; # Copy deformed endocardial surfaces (interface nodes) from wall region to cavity region

# Update locations of cavity nodes based on average positions of endocardial 'rings' (in 1 represents the RC coordinate along which node 118 is averaged among nodes 14..17)
fem update solution cavity_reference average 118 in 1 node 14..17 region $LV_CAVITY;

####################################### Calculates LV Cavity Volume #########################################
# Output current pressure
read command;Extract_Pressure;

# List the deformed cavity volume
fem list elements;output_cavity_volume/LVCavity_${TCa} deformed total region $LV_CAVITY;
$LV_CAVITY_VOLUME=`awk '\$1=="Total"{printf("%f",\$4);exit}\' output_cavity_volume/LVCavity_${TCa}.opelem`;
print "\nLV cavity volume   (Current) = ${LV_CAVITY_VOLUME}\n";
print "  LV cavity pressure (Current) = ${LV_P} kPa\n";
print "  LV Activation      (Current) = ${TCa} kPa \n";

#############################################################################################################
########## Redefine files required for the wall mechanics simulation ############################
fem define coupling;r;coupled_reset
