gfx create material heart ambient 0.3 0 0.3 diffuse 1 0 0 specular 0.5 0.5 0.5 shininess 0.5;
gfx create material bluey ambient 0 0.25 0.5 diffuse 0 0.4 1 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3
gfx create material copper ambient 1 0.2 0 diffuse 0.6 0.3 0 emission 0 0 0 specular 0.7 0.7 0.5 alpha 1 shininess 0.3
gfx create material gold ambient 1 0.4 0 diffuse 1 0.7 0 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.8
gfx create material silver ambient 0.4 0.4 0.4 diffuse 0.7 0.7 0.7 emission 0 0 0 specular 0.7 0.7 0.7 alpha 1 shininess 0.6
gfx cre mat trans_purple ambient 0.4 0 0.9 diffuse 0.4 0 0.9 alpha 0.3
gfx create spectrum error_spectrum;
gfx modify spectrum error_spectrum clear overwrite_colour;
gfx modify spectrum error_spectrum linear range 0.0007865 3 extend_above extend_below white_to_red colour_range 0 1 component 1;

gfx cre win 1;
gfx mod win 1 back colour "0,0,0";
gfx create material heart ambient 0.3 0 0.3 diffuse 1 0 0 specular 0.5 0.5 0.5 shininess 0.5;

$DS = 23;
$tot = 30;
$ED_frame=1;

$time=0;

for ($i=$DS+1;$i<=$tot;$i=$i+1) 
{
	$file="ForwardSolveExfile/LVInflation_".${i};	
	gfx read node;$file region LVModel time $time;
	$time++;	
}

$i=$ED_frame;
$file="ForwardSolveExfile/LVInflation_".${i};	
gfx read node;$file region LVModel time $time;
gfx read element;"ForwardSolveExfile/LVInflation_".${i} region LVModel
gfx  modify g_element LVModel general clear circle_discretization 48 default_coordinate deformed element_discretization "12*12*12" native_discretization none
gfx mod g_e LVModel surfaces exterior face xi3_0 mat muscle
gfx mod g_elem LVModel cylinders constant_radius 0.2 material muscle;

$time=0;
for ($i=$DS+1;$i<=$tot;$i=$i+1)  
{
	$file="ForwardSolveExfile/Surface_epi_".${i};
	gfx read data;$file region Surface_Points_Epi time $time;
	$time++;
	
}
$i=$ED_frame;
$file="ForwardSolveExfile/Surface_epi_".${i};
gfx read data;$file region Surface_Points_Epi time $time;
gfx modify g_element Surface_Points_Epi data_points glyph sphere general size "2*2*2" centre 0,0,0 font default select_on material green selected_material default_selected;

$time=0;
for ($i=$DS+1;$i<=$tot;$i=$i+1)  
{
	$file="ForwardSolveError/EpiError_".${i};
	gfx read data;$file region CAPEDEpi_Error_ToED time $time;
	$time++;
}
$i=$ED_frame;
$file="ForwardSolveError/EpiError_".${i};
gfx read data;$file region CAPEDEpi_Error_ToED time $time;
gfx modify g_element CAPEDEpi_Error_ToED general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element CAPEDEpi_Error_ToED lines select_on material gold selected_material default_selected;
gfx modify g_element CAPEDEpi_Error_ToED data_points glyph arrow_solid general size "0.8*0.8*0.8" centre 0,0,0 font default orientation error scale_factors "1*0.5*0.5" select_on material gold selected_material default_selected;

$time=0;
for ($i=$DS+1;$i<=$tot;$i=$i+1) 
{
	$file="ForwardSolveExfile/Surface_endo_".${i};	
	gfx read data;$file region Surface_Points_Endo time $time;
	$time++;
	
}

$i=$ED_frame;
$file="Surface_endo_".${i};	
gfx read data;$file region Surface_Points_Endo time $time;
gfx modify g_element Surface_Points_Endo data_points glyph sphere general size "2*2*2" centre 0,0,0 font default select_on material gold selected_material default_selected;

$time=0;
for ($i=$DS+1;$i<=$tot;$i=$i+1)  
{
	$file="ForwardSolveError/EndoError_".${i};
	gfx read data;$file region CAPEDEndo_Error_ToED time $time;
	$time++;
}

$i=$ED_frame;
$file="ForwardSolveError/EndoError_".${i};
gfx read data;$file region CAPEDEndo_Error_ToED time $time;

gfx modify g_element CAPEDEndo_Error_ToED general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element CAPEDEndo_Error_ToED lines select_on material gold selected_material default_selected;
gfx modify g_element CAPEDEndo_Error_ToED data_points glyph arrow_solid general size "0.8*0.8*0.8" centre 0,0,0 font default orientation error scale_factors "1*0.5*0.5" select_on material gold selected_material default_selected;


gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout front_back ortho_axes -x -z eye_spacing 0.25 width 1475 height 855;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 23.6889 11.0331 -167.865 interest_point 23.6889 11.0331 9.3848 up_vector -1 0 -0 view_angle 57.3856 near_clipping_plane 1.7725 far_clipping_plane 633.431 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set current_pane 2;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 23.6889 11.0331 186.635 interest_point 23.6889 11.0331 9.3848 up_vector -1 0 0 view_angle 57.3856 near_clipping_plane 2.34198 far_clipping_plane 836.945 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

