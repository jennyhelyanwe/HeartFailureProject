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

gfx read node;output/LV_EndSystole_OptTCa
gfx read elem;output/LV_EndSystole_OptTCa
gfx modify g_element LVEndSystole_OptTCa general clear circle_discretization 6 default_coordinate deformed element_discretization "12*12*12" native_discretization none;
gfx modify g_element LVEndSystole_OptTCa lines select_on material default selected_material default_selected;
gfx modify g_element LVEndSystole_OptTCa cylinders constant_radius 0.3 select_on material muscle selected_material default_selected render_shaded;
gfx modify g_element LVEndSystole_OptTCa surfaces face xi3_1 select_on material muscle selected_material default_selected render_shaded;

gfx read data;CAPES_Endo
gfx modify g_element CAPES_Endo general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element CAPES_Endo lines select_on material default selected_material default_selected;
gfx modify g_element CAPES_Endo data_points glyph sphere general size "1*1*1" centre 0,0,0 font default select_on material green selected_material default_selected;

gfx change data_offset 10000
gfx read data;CAPES_Epi
gfx modify g_element CAPES_Epi general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element CAPES_Epi lines select_on material default selected_material default_selected;
gfx modify g_element CAPES_Epi data_points glyph sphere general size "1*1*1" centre 0,0,0 font default select_on material green selected_material default_selected;

gfx change data_offset 10000
gfx read data;CAPESEndo_Error_ToES
gfx modify g_element CAPESEndo_Error_ToES general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element CAPESEndo_Error_ToES lines select_on material gold selected_material default_selected;
gfx modify g_element CAPESEndo_Error_ToES data_points glyph arrow_solid general size "0.8*0.8*0.8" centre 0,0,0 font default orientation error scale_factors "1*0.5*0.5" select_on material gold selected_material default_selected;

gfx change data_offset 10000
gfx read data;CAPESEpi_Error_ToES
gfx modify g_element CAPESEpi_Error_ToES general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element CAPESEpi_Error_ToES lines select_on material gold selected_material default_selected;
gfx modify g_element CAPESEpi_Error_ToES data_points glyph arrow_solid general size "0.8*0.8*0.8" centre 0,0,0 font default orientation error scale_factors "1*0.5*0.5" select_on material gold selected_material default_selected;

gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 1475 height 973;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point -35.4206 -168.401 -199.274 interest_point 16.3439 10.0188 3.15872 up_vector -0.980989 0.0888717 0.17252 view_angle 31.99 near_clipping_plane 2.74759 far_clipping_plane 981.893 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

