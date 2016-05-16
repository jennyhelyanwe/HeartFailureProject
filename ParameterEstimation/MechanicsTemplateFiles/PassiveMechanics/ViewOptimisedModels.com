gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear range 0 0.15 extend_above white_to_red colour_range 0 1 component 1;
gfx modify spectrum default linear reverse range -0.15 0 extend_below white_to_blue colour_range 0 1 component 1;
gfx create material black normal_mode ambient 0 0 0 diffuse 0 0 0 emission 0 0 0 specular 0.3 0.3 0.3 alpha 1 shininess 0.2;
gfx create material blue normal_mode ambient 0 0 0.5 diffuse 0 0 1 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.2;
gfx create material bone normal_mode ambient 0.7 0.7 0.6 diffuse 0.9 0.9 0.7 emission 0 0 0 specular 0.1 0.1 0.1 alpha 1 shininess 0.2;
gfx create material default normal_mode ambient 1 1 1 diffuse 1 1 1 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx create material default_selected normal_mode ambient 1 0.2 0 diffuse 1 0.2 0 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx create material gold normal_mode ambient 1 0.4 0 diffuse 1 0.7 0 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3;
gfx create material gray50 normal_mode ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5 emission 0.5 0.5 0.5 specular 0.5 0.5 0.5 alpha 1 shininess 0.2;
gfx create material green normal_mode ambient 0 0.5 0 diffuse 0 1 0 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.1;
gfx create material muscle normal_mode ambient 0.4 0.14 0.11 diffuse 0.5 0.12 0.1 emission 0 0 0 specular 0.3 0.5 0.5 alpha 1 shininess 0.2;
gfx create material red normal_mode ambient 0.5 0 0 diffuse 1 0 0 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.2;
gfx create material silver normal_mode ambient 0.4 0.4 0.4 diffuse 0.7 0.7 0.7 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3;
gfx create material tissue normal_mode ambient 0.9 0.7 0.5 diffuse 0.9 0.7 0.5 emission 0 0 0 specular 0.2 0.2 0.3 alpha 1 shininess 0.2;
gfx create material transparent_gray50 normal_mode ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5 emission 0.5 0.5 0.5 specular 0.5 0.5 0.5 alpha 0 shininess 0.2;
gfx create material white normal_mode ambient 1 1 1 diffuse 1 1 1 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx create spectrum error_spectrum;
gfx modify spectrum error_spectrum clear overwrite_colour;
gfx modify spectrum error_spectrum linear range 0.0007865 3 extend_above extend_below white_to_red colour_range 0 1 component 1;


######################## Visualise inflated models ##################
$DS = 23
$tot = 30
$ED = 1

$t=0;
for ($i=$DS+1;$i<=$tot;$i=$i+1)
{
	gfx read node;OptimisedExfile/LVInflation_$i time $t;
	$t ++;
}
gfx read node;OptimisedExfile/LVInflation_$ED time $t;
gfx read element;OptimisedExfile/LVInflation_$ED;

gfx modify g_element LVInflation_Update general clear circle_discretization 6 default_coordinate deformed element_discretization "12*12*12" native_discretization none;
gfx modify g_element LVInflation_Update lines select_on material default selected_material default_selected;
gfx modify g_element LVInflation_Update cylinders constant_radius 0.3 select_on material muscle selected_material default_selected render_shaded;
gfx modify g_element LVInflation_Update surfaces face xi3_0 select_on material muscle selected_material default_selected render_shaded;
gfx modify g_element "/" points domain_point tessellation default_points LOCAL glyph axes_solid_colour size "20*20*20" offset 0,0,0 font default select_on material default selected_material default_selected render_shaded;

$t=0;
for ($i=$DS+1;$i<=$tot;$i=$i+1)
{
	gfx read data;OptimisedExfile/Surface_epi_$i time $t region Surface_Points_Epi ;
	$t ++;
	
}
gfx read data;OptimisedExfile/Surface_epi_$ED  time $t region Surface_Points_Epi;
gfx modify g_element Surface_Points_Epi general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element Surface_Points_Epi data_points glyph sphere general size "1*1*1" centre 0,0,0 font default select_on material blue selected_material default_selected;

$t=0;
for ($i=$DS+1;$i<=$tot;$i=$i+1)
{
	gfx read data;OptimisedExfile/Surface_endo_$i time $t region Surface_Points_Endo;
	$t ++;
	
}
gfx read data;OptimisedExfile/Surface_endo_$ED   time $t region Surface_Points_Endo;
gfx modify g_element Surface_Points_Endo general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element Surface_Points_Endo data_points glyph sphere general size "1*1*1" centre 0,0,0 font default select_on material blue selected_material default_selected;

################### Visualising error projections ##################
gfx change data_offset 10000
$t=0;
for ($i=$DS+1;$i<=$tot;$i=$i+1)
{
	gfx read data;OptimisedError/EndoError_$i time $t region EndoError;
	$t ++;
	
}
gfx read data;OptimisedError/EndoError_$ED time $t region EndoError;
gfx def field EndoError/error_mag magnitude field error;
gfx modify g_element /EndoError/ general clear;
gfx modify g_element /EndoError/ lines domain_mesh1d coordinate coordinates tessellation temp4 LOCAL line line_base_size 0 select_on material default data error spectrum default selected_material default_selected render_shaded;
gfx modify g_element /EndoError/ points domain_datapoints coordinate coordinates tessellation default_points LOCAL glyph arrow_solid size "0.8*0.8*0.8" offset 0,0,0 font default orientation error scale_factors "1*0.5*0.5" select_on material gold data error_mag spectrum error_spectrum selected_material default_selected render_shaded;

$t=0;
for ($i=$DS+1;$i<=$tot;$i=$i+1)
{
	gfx read data;OptimisedError/EpiError_$i time $t region EpiError;
	$t ++;
	
}
gfx read data;OptimisedError/EpiError_$ED time $t region EpiError;
gfx def field EpiError/error_mag magnitude field error;
gfx modify g_element /EpiError/ general clear;
gfx modify g_element /EpiError/ lines domain_mesh1d coordinate coordinates tessellation temp4 LOCAL line line_base_size 0 select_on material gold selected_material default_selected render_shaded;
gfx modify g_element /EpiError/ points domain_datapoints coordinate coordinates tessellation default_points LOCAL glyph arrow_solid size "0.8*0.8*0.8" offset 0,0,0 font default orientation error scale_factors "1*0.5*0.5" select_on material gold data error_mag spectrum error_spectrum selected_material default_selected render_shaded;




gfx create window 1 double_buffer;
gfx modify window 1 image scene "/" filter default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 1475 height 855;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point -51.6351 -224.288 -262.683 interest_point 16.3439 10.0188 3.15872 up_vector -0.980988 0.0888722 0.172521 view_angle 31.99 near_clipping_plane 88.8116 far_clipping_plane 1067.96 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

$t = 0;
for ($i=$DS+1;$i<=$tot;$i=$i+1)
{
	gfx print jpg file OptimisedExfile/Image_$i window 1;
	$t ++;
	gfx set time $t;	
}
gfx set time $t;
gfx print jpg file OptimisedExfile/Image_$ED window 1;

####### Calculate strain
gfx read node UPFInitialModelWithFibre
gfx read elem UPFInitialModelWithFibre

gfx define field F gradient coordinate coordinates field deformed
gfx define field F_transpose transpose source_number_of_rows 3 field F
gfx define field identity3 composite 1 0 0 0 1 0 0 0 1

gfx define field C matrix_multiply number_of_rows 3 fields F_transpose F
gfx define field E2 add fields C identity3 scale_factors 1 -1
gfx define field E coordinate_system rectangular_cartesian scale field E2 scale_factors 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5

gfx define field principal_strains eigenvalues field E
gfx define field principal_strain_vectors eigenvectors eigenvalues principal_strains
gfx define field deformed_principal_strain_vectors matrix_multiply number_of_rows 3 fields principal_strain_vectors F_transpose

gfx define field deformed_principal_strain_vector1 composite deformed_principal_strain_vectors.1 deformed_principal_strain_vectors.2 deformed_principal_strain_vectors.3
gfx define field deformed_principal_strain_vector2 composite deformed_principal_strain_vectors.4 deformed_principal_strain_vectors.5 deformed_principal_strain_vectors.6
gfx define field deformed_principal_strain_vector3 composite deformed_principal_strain_vectors.7 deformed_principal_strain_vectors.8 deformed_principal_strain_vectors.9
# since above vectors have the stretch as their magnitude, normalize them:
gfx define field norm_def_principal_strain_vector1 normalise field deformed_principal_strain_vector1
gfx define field norm_def_principal_strain_vector2 normalise field deformed_principal_strain_vector2
gfx define field norm_def_principal_strain_vector3 normalise field deformed_principal_strain_vector3

gfx define field principal_strain1 composite principal_strains.1
gfx define field principal_strain2 composite principal_strains.2
gfx define field principal_strain3 composite principal_strains.3

#Calculate the deformed fibre axes
gfx define field fibre_axes fibre_axes coordinate coordinates fibre fibres
gfx define field deformed_fibre_axes matrix_multiply number_of_rows 3 fields fibre_axes F_transpose

gfx define field deformed_fibre composite deformed_fibre_axes.1 deformed_fibre_axes.2 deformed_fibre_axes.3
gfx define field deformed_sheet composite deformed_fibre_axes.4 deformed_fibre_axes.5 deformed_fibre_axes.6
gfx define field norm_def_fibre normalise field deformed_fibre
gfx define field def_fibre_cross cross_product dimension 3 fields norm_def_fibre deformed_sheet
gfx define field norm_def_fibre_cross normalise field def_fibre_cross
gfx define field norm_def_fibre_cross_normal cross_product dimension 3 fields norm_def_fibre norm_def_fibre_cross
gfx define field orthonormal_deformed_fibre_axes composite norm_def_fibre norm_def_fibre_cross_normal norm_def_fibre_cross

#gfx modify g_element "/" points domain_mesh_highest_dimension coordinate deformed tessellation default_points LOCAL glyph sphere size "5*0.5*0.5" offset 0,0,0 font default orientation norm_def_principal_strain_vector1 scale_factors "1*1*1" cell_centres select_on material default data principal_strain1 spectrum default selected_material default_selected render_shaded;

#gfx define tessellation default_points minimum_divisions "6" refinement_factors "1" circle_divisions 12;
#gfx define tessellation temp3 minimum_divisions "12*12*12" refinement_factors "1" circle_divisions 12;
#gfx define tessellation temp4 minimum_divisions "4*4*4" refinement_factors "1" circle_divisions 12;
