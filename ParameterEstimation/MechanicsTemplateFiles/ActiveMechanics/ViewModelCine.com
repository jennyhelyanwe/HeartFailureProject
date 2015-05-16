gfx create material heart ambient 0.3 0 0.3 diffuse 1 0 0 specular 0.5 0.5 0.5 shininess 0.5;
gfx create material bluey ambient 0 0.25 0.5 diffuse 0 0.4 1 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3
gfx create material copper ambient 1 0.2 0 diffuse 0.6 0.3 0 emission 0 0 0 specular 0.7 0.7 0.5 alpha 1 shininess 0.3
gfx create material gold ambient 1 0.4 0 diffuse 1 0.7 0 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.8
gfx create material silver ambient 0.4 0.4 0.4 diffuse 0.7 0.7 0.7 emission 0 0 0 specular 0.7 0.7 0.7 alpha 1 shininess 0.6
gfx cre mat trans_purple ambient 0.4 0 0.9 diffuse 0.4 0 0.9 alpha 0.3

gfx cre win 1;
gfx mod win 1 back colour "0,0,0";
gfx create material heart ambient 0.3 0 0.3 diffuse 1 0 0 specular 0.5 0.5 0.5 shininess 0.5;


for ($i=1;$i<=30;$i=$i+1) 
{
	gfx read node;outputAll/LVModel_$i time $i;
	
}

gfx read element;outputAll/LVModel_1
gfx  modify g_element LVModel general clear circle_discretization 48 default_coordinate deformed element_discretization "12*12*12" native_discretization none
gfx mod g_e LVModel surfaces exterior face xi3_0 mat muscle

gfx mod g_elem LVModel cylinders constant_radius 0.2 material muscle;

gfx read data;CAPES_Endo
gfx modify g_element CAPES_Endo general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element CAPES_Endo lines select_on material default selected_material default_selected;
gfx modify g_element CAPES_Endo data_points glyph sphere general size "2*2*2" centre 0,0,0 font default select_on material green selected_material default_selected;

gfx change data_offset 10000
gfx read data;CAPES_Epi
gfx modify g_element CAPES_Epi general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element CAPES_Epi lines select_on material default selected_material default_selected;
gfx modify g_element CAPES_Epi data_points glyph sphere general size "2*2*2" centre 0,0,0 font default select_on material green selected_material default_selected;

gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout front_back ortho_axes -x -z eye_spacing 0.25 width 1475 height 857;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 23.6889 11.0331 -167.865 interest_point 23.6889 11.0331 9.3848 up_vector -1 0 -0 view_angle 52.7152 near_clipping_plane 1.7725 far_clipping_plane 633.431 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set current_pane 2;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 23.6889 11.0331 186.635 interest_point 23.6889 11.0331 9.3848 up_vector -1 0 0 view_angle 52.7152 near_clipping_plane 2.34198 far_clipping_plane 836.945 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;
