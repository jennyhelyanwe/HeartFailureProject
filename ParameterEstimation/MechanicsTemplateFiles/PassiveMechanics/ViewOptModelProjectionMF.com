gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0 1 extend_above extend_below rainbow colour_range 0 1 component 1;
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

######################## Visualise inflated models ##################
$t = 0;
for ($i=24;$i<=30;$i=$i+1)
{
	gfx read node;OptimisedExfile/LVInflation_$i time $t;
	$t ++;
	
}
gfx read node;OptimisedExfile/LVInflation_1 time $t;
gfx read element;OptimisedExfile/LVInflation_1;

gfx modify g_element LVInflation_Update general clear circle_discretization 6 default_coordinate deformed element_discretization "12*12*12" native_discretization none;
gfx modify g_element LVInflation_Update lines select_on material default selected_material default_selected;
gfx modify g_element LVInflation_Update cylinders constant_radius 0.3 select_on material muscle selected_material default_selected render_shaded;
gfx modify g_element LVInflation_Update surfaces face xi3_1 select_on material muscle selected_material default_selected render_shaded;

##################### Visualised ground truth surface points ############

$t = 0;
for ($i=24;$i<=30;$i=$i+1)
{
	gfx read data;OptimisedExfile/Surface_endo_reg_$i time $t region endoReg;
	$t ++;
	
}
gfx read data;OptimisedExfile/Surface_endo_reg_1 time $t region endoReg;
gfx modify g_element endoReg/trans general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element endoReg/trans data_points glyph sphere general size "1*1*1" centre 0,0,0 font default select_on material green selected_material default_selected;

$t = 0;
for ($i=24;$i<=30;$i=$i+1) 
{
	gfx read data;OptimisedExfile/Surface_epi_reg_$i time $t region epiReg;
	$t ++;
	
}
gfx read data;OptimisedExfile/Surface_epi_reg_1 time $t region epiReg;
gfx modify g_element epiReg/trans general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element epiReg/trans data_points glyph sphere general size "1*1*1" centre 0,0,0 font default select_on material green selected_material default_selected;

$t = 0;
for ($i=24;$i<=30;$i=$i+1)
{
	gfx read data;OptimisedExfile/Surface_epi_nonreg_$i time $t region Surface_Points_Epi ;
	$t ++;
	
}
gfx read data;OptimisedExfile/Surface_epi_nonreg_1  time $t region Surface_Points_Epi;
gfx modify g_element Surface_Points_Epi general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element Surface_Points_Epi data_points glyph sphere general size "1*1*1" centre 0,0,0 font default select_on material blue selected_material default_selected;

$t = 0;
for ($i=24;$i<=30;$i=$i+1)
{
	gfx read data;OptimisedExfile/Surface_endo_nonreg_$i time $t region Surface_Points_Endo;
	$t ++;
	
}
gfx read data;OptimisedExfile/Surface_endo_nonreg_1   time $t region Surface_Points_Endo;
gfx modify g_element Surface_Points_Endo general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element Surface_Points_Endo data_points glyph sphere general size "1*1*1" centre 0,0,0 font default select_on material blue selected_material default_selected;

################### Visualising error projections ##################
gfx change data_offset 10000
$t = 0;
for ($i=24;$i<=30;$i=$i+1)
{
	gfx read data;OptimisedError/EndoError_$i time $t region EndoError;
	$t ++;
	
}
gfx read data;OptimisedError/EndoError_1 time $t region EndoError;

gfx modify g_element EndoError general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element EndoError lines select_on material gold selected_material default_selected;
gfx modify g_element EndoError data_points glyph arrow_solid general size "0.8*0.8*0.8" centre 0,0,0 font default orientation error scale_factors "1*0.5*0.5" select_on material gold selected_material default_selected;

$t = 0;
for ($i=24;$i<=30;$i=$i+1)
{
	gfx read data;OptimisedError/EpiError_$i time $t region EpiError;
	$t ++;
	
}
gfx read data;OptimisedError/EpiError_1 time $t region EpiError;
gfx modify g_element EpiError general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element EpiError lines select_on material gold selected_material default_selected;
gfx modify g_element EpiError data_points glyph arrow_solid general size "0.8*0.8*0.8" centre 0,0,0 font default orientation error scale_factors "1*0.5*0.5" select_on material gold selected_material default_selected;



gfx create window 1 double_buffer;
gfx modify window 1 image scene "/" filter default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 1475 height 855;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point -51.6351 -224.288 -262.683 interest_point 16.3439 10.0188 3.15872 up_vector -0.980988 0.0888722 0.172521 view_angle 31.99 near_clipping_plane 88.8116 far_clipping_plane 1067.96 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

$t = 0;
for ($i=24;$i<=30;$i=$i+1)
{
	gfx print jpg file OptimisedExfile/Image_$i window 1;
	$t ++;
	gfx set time $t;	
}
gfx set time $t;
gfx print jpg file OptimisedExfile/Image_1 window 1;


