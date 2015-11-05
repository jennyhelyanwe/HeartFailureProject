
$DS = 16;
$ED = 1;
$tot = 25;
$study = STF_01
$t = 0
for ($i = $ED;$i<=$DS; $i++){
    gfx read data ${study}/Active/${study}_Surface_Points_Epi_$i time $t region epi;
    $t = $t + 1;
}



$t = 0
for ($i = $ED;$i<=$DS; $i++){
    gfx read data ${study}/Active/${study}_Surface_Points_Endo_$i time $t region endo;
    $t = $t + 1;
}




$t = 0
for ($i = $ED;$i<=$DS; $i++){
    gfx read node ${study}_$i time $t region model;
    $t = $t + 1;
}
gfx read elem ${study}_$DS region model;

gfx modify g_element "/" general clear;
gfx modify g_element /epi/ general clear;
gfx modify g_element /epi/ data_points coordinate coordinates LOCAL glyph sphere_hires general size "1*1*1" centre 0,0,0 font default select_on material green selected_material default_selected;
gfx modify g_element /endo/ general clear;
gfx modify g_element /endo/ data_points coordinate coordinates LOCAL glyph sphere_hires general size "1*1*1" centre 0,0,0 font default select_on material gold selected_material default_selected;
gfx modify g_element /model/ general clear;
gfx modify g_element /model/ cylinders coordinate coordinates tessellation default LOCAL circle_discretization 6 constant_radius 0.2 select_on material green selected_material default_selected render_shaded;
gfx modify g_element /model/ surfaces coordinate coordinates face xi3_0 tessellation default LOCAL select_on material green selected_material default_selected render_shaded;

gfx define tessellation default minimum_divisions "10" refinement_factors "4";

gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 635 height 855;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 19.0869 17.5518 -243.894 interest_point 10.787 -2.05385 -2.27965 up_vector -0.998751 0.0390837 -0.0311374 view_angle 40 near_clipping_plane 2.42551 far_clipping_plane 866.793 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;


