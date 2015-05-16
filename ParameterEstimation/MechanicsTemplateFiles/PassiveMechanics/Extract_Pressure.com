################################################# Extract Pressure #########################################################
fem list variable;pressure/pressure # Creates an .opvari file containing values for nodal and auxiliary degrees of freedom
$LINE=1007; # Specifies a line which represents pressure acting on the endocardial face of an element
$LV_P=`sed -e s%D%E% pressure/pressure.opvari | awk -v line=$LINE 'NR==line{printf("%.2f",\$4)}'` # Extracts pressure value

print "\033[0;30;42m    Current pressure = $LV_P kPa    \033[0m\n";

############################################################################################################################
