gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0 100 extend_above extend_below rainbow colour_range 0 1 component 1;

gfx read node OptimisedExfile/CurrentContracted_12 region model
gfx read elem OptimisedExfile/CurrentContracted_12 region model
gfx read data TotalStress2D_12

gfx modify g_element "/" general clear;
gfx modify g_element "/" points domain_datapoints coordinate coordinates tessellation default_points LOCAL glyph sphere size "4*4*4" offset 0,0,0 font default select_on material default data total_mps spectrum default selected_material default_selected render_shaded;
gfx modify g_element /model/ general clear;
gfx modify g_element /model/ lines domain_mesh1d coordinate deformed tessellation default LOCAL circle_extrusion line_base_size 0.3 select_on material muscle selected_material default_selected render_shaded;
#gfx modify g_element /model/ points domain_point tessellation default_points LOCAL glyph axes_solid_colour size "20*20*20" offset 0,0,0 font default select_on material default selected_material default_selected render_shaded;
#gfx modify g_element /model/ surfaces domain_mesh2d coordinate deformed face xi2_0 tessellation default LOCAL select_on material muscle selected_material default_selected render_shaded;

gfx define tessellation default minimum_divisions "10" refinement_factors "4" circle_divisions 12;
gfx define tessellation default_points minimum_divisions "1" refinement_factors "1" circle_divisions 12;

gfx destroy lines 1..29,33,36,39,41..69,72,74,76..81 group model
#gfx destroy faces 3,7,11,14,19,23,27,30 group model

gfx create window 1 double_buffer;
gfx modify window 1 image scene "/" filter default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout 2d ortho_axes x z eye_spacing 0.25 width 635 height 974;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 2.32434e-310 6.95322e-310 3.08043e-316 texture none;
gfx modify window 1 view parallel eye_point 194.644 0.0479511 0.943104 interest_point 2.82773 0.0479511 0.943104 up_vector -0 -0 -1 view_angle 40 near_clipping_plane 0.183408 far_clipping_plane 836.974 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

gfx print png file TotalStress2D.png window 1;
q
