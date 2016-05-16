gfx create spectrum error_spectrum;
gfx modify spectrum error_spectrum clear overwrite_colour;
gfx modify spectrum error_spectrum linear range 0.0007865 3 extend_above extend_below white_to_red colour_range 0 1 component 1;

$DS = 23
$ED = 1
$tot = 30

gfx read node DSWall
gfx read elem DSWall

$t = 0
for ($i = $DS; $i <=$tot; $i=$i+1){
    gfx read data DisplacementMSE/EpiData_$i time $t region data_epi;
    $t = $t + 1
}
gfx read data DisplacementMSE/EpiData_$ED time $t region data_epi;

$t = 0
for ($i = $DS; $i <=$tot; $i=$i+1){
    gfx read data DisplacementMSE/EpiDisp_$i time $t region disp_epi;
    $t = $t + 1
}
gfx read data DisplacementMSE/EpiDisp_$ED time $t region disp_epi;


gfx modify g_element "/" general clear;
gfx modify g_element "/" surfaces coordinate coordinates face xi3_0 tessellation default LOCAL select_on material muscle selected_material default_selected render_shaded;
gfx modify g_element "/" lines coordinate coordinates tessellation default LOCAL select_on material default selected_material default_selected;
gfx modify g_element "/" cylinders constant_radius 0.3 select_on material muscle selected_material default_selected render_shaded;
gfx modify g_element /data_epi/ general clear;
gfx modify g_element /data_epi/ data_points coordinate coordinates LOCAL glyph sphere general size "1*1*1" centre 0,0,0 font default select_on material green selected_material default_selected;
gfx def field /disp_epi/error_mag magnitude field error;
gfx modify g_element /disp_epi/ general clear;
gfx modify g_element /disp_epi/ points domain_datapoints coordinate coordinates tessellation default_points LOCAL glyph arrow_solid size "0.8*0.8*0.8" offset 0,0,0 font default orientation error scale_factors "1*0.5*0.5" select_on material gold data error_mag spectrum error_spectrum selected_material default_selected render_shaded;

$t = 0
for ($i = $DS; $i <=$tot; $i=$i+1){
    gfx read data DisplacementMSE/EndoData_$i time $t region data_endo;
    $t = $t + 1
}
gfx read data DisplacementMSE/EndoData_$ED time $t region data_endo;

$t = 0
for ($i = $DS; $i <=$tot; $i=$i+1){
    gfx read data DisplacementMSE/EndoDisp_$i time $t region disp_endo;
    $t = $t + 1
}
gfx read data DisplacementMSE/EndoDisp_$ED time $t region disp_endo;

gfx modify g_element /data_endo/ general clear;
gfx modify g_element /data_endo/ data_points coordinate coordinates LOCAL glyph sphere general size "1*1*1" centre 0,0,0 font default select_on material green selected_material default_selected;
gfx def field /disp_endo/error_mag magnitude field error;
gfx modify g_element /disp_endo/ general clear;
gfx modify g_element /disp_endo/ points domain_datapoints coordinate coordinates tessellation default_points LOCAL glyph arrow_solid size "0.8*0.8*0.8" offset 0,0,0 font default orientation error scale_factors "1*0.5*0.5" select_on material gold data error_mag spectrum error_spectrum selected_material default_selected render_shaded;

gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 512 height 512;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 15.6842 34.3586 -319.929 interest_point 18.3685 -2.666 1.3085 up_vector -0.99326 0.113913 0.0214291 view_angle 40 near_clipping_plane 3.23375 far_clipping_plane 1155.63 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;
