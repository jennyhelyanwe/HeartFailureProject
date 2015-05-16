gfx create material heart ambient 0.3 0 0.3 diffuse 1 0 0 specular 0.5 0.5 0.5 shininess 0.5;
gfx create material bluey ambient 0 0.25 0.5 diffuse 0 0.4 1 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3
gfx create material copper ambient 1 0.2 0 diffuse 0.6 0.3 0 emission 0 0 0 specular 0.7 0.7 0.5 alpha 1 shininess 0.3
gfx create material gold ambient 1 0.4 0 diffuse 1 0.7 0 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.8
gfx create material silver ambient 0.4 0.4 0.4 diffuse 0.7 0.7 0.7 emission 0 0 0 specular 0.7 0.7 0.7 alpha 1 shininess 0.6
gfx cre mat trans_purple ambient 0.4 0 0.9 diffuse 0.4 0 0.9 alpha 0.3

gfx cre win 1;
gfx mod win 1 back colour "0,0,0";
gfx create material heart ambient 0.3 0 0.3 diffuse 1 0 0 specular 0.5 0.5 0.5 shininess 0.5;


for ($i=0;$i<=10;$i=$i+1) 
{
	gfx read node;output/LV_Inflation_$i time $i;
	
}

gfx read element;output/LV_Inflation_0

gfx modify g_element LVInflation general clear circle_discretization 48 default_coordinate deformed element_discretization "12*12*12" native_discretization none;
gfx modify g_element LVInflation cylinders constant_radius 0.2 select_on material muscle selected_material default_selected render_shaded;
gfx modify g_element LVInflation surfaces face xi3_0 select_on material muscle selected_material default_selected render_shaded;

gfx read data;CAPED_Endo
gfx modify g_element CAPED_Endo general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element CAPED_Endo lines select_on material default selected_material default_selected;
gfx modify g_element CAPED_Endo data_points glyph sphere general size "2*2*2" centre 0,0,0 font default select_on material green selected_material default_selected;

gfx change data_offset 10000
gfx read data;CAPED_Epi
gfx modify g_element CAPED_Epi general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element CAPED_Epi lines select_on material default selected_material default_selected;
gfx modify g_element CAPED_Epi data_points glyph sphere general size "2*2*2" centre 0,0,0 font default select_on material green selected_material default_selected;


gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout front_back ortho_axes -x -z eye_spacing 0.25 width 1475 height 857;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 23.6889 11.0331 -167.865 interest_point 23.6889 11.0331 9.3848 up_vector -1 0 0 view_angle 40.3678 near_clipping_plane 1.7725 far_clipping_plane 633.431 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set current_pane 2;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 23.6889 11.0331 186.635 interest_point 23.6889 11.0331 9.3848 up_vector -1 0 0 view_angle 40.3678 near_clipping_plane 2.34198 far_clipping_plane 836.945 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

