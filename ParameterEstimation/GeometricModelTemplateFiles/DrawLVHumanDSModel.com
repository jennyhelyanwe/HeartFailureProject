# Create one window to display the fitted mesh to -axis data
gfx cre window 1;


# Read in the node file
gfx read nodes DSWall

# Read in the element file
gfx read elements DSWall

gfx  modify g_element DSWall general clear circle_discretization 48 default_coordinate coordinates element_discretization "12*12*12" native_discretization none;

gfx modify g_element DSWall node_points glyph sphere general size "2*2*2" centre 0,0,0 material default;
gfx modify g_element DSWall element_points; 
gfx modify g_element DSWall cylinders constant_radius 0.2 material muscle;
gfx modify g_element DSWall surfaces material muscle;


gfx change node_offset 1000
gfx change face_offset 1000
gfx change element_offset 1000


## Read in the node file
gfx read nodes DSCavity

## Read in the element file
gfx read elements DSCavity

gfx  modify g_element DSCavity general clear circle_discretization 48 default_coordinate coordinates element_discretization "12*12*12" native_discretization none;

gfx modify g_element DSCavity node_points glyph sphere general size "2*2*2" centre 0,0,0 material default;
gfx modify g_element DSCavity element_points; 
gfx modify g_element DSCavity cylinders constant_radius 0.2 material muscle;
gfx modify g_element DSCavity surfaces material muscle;

# Change some visualization paramaters
gfx modify window 1 image view_all;
gfx mod win 1 background colour 1 1 1
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout 2d ortho_axes -z x eye_spacing 0.25 width 1473 height 973;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 19.3275 0.894542 -258.998 interest_point 19.3275 0.894542 16.5363 up_vector -1 0 0 view_angle 41.8673 near_clipping_plane 2.75534 far_clipping_plane 984.664 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;


