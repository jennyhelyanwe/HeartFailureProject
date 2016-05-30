gfx create material heart ambient 0.3 0 0.3 diffuse 1 0 0 specular 0.5 0.5 0.5 shininess 0.5;
gfx cre mat trans_purple ambient 0.4 0 0.9 diffuse 0.4 0 0.9 alpha 0.3;

gfx cre win;
gfx mod win 1 background colour 1 1 1;
gfx edit scene;

# Set up the output directory



for ($i=1;$i<=20;$i=$i+1)
{	
	
	$FILENAME="model/SCD0001301_$i.model.exnode";	
	gfx read node ${FILENAME} time $i;
	gfx mod g_elem heart node_points glyph sphere size 2;
}


gfx read element GlobalHermiteParam;

gfx modify g_element heart node_points glyph sphere general size "2*2*2" centre 0,0,0 select_on material default;
gfx  Mod g_element heart general clear circle_discretization 6 default_coordinate coordinates element_discretization "12*12*12" native_discretization none;

gfx modify g_element heart cylinders constant_radius 0.3;
gfx  Mod g_e heart surfaces exterior face xi3_0 mat heart
gfx  Mod g_e heart surfaces exterior face xi3_1 mat trans_purple

gfx modify window 1 image view_all;

gfx modify window 1 layout 2d ortho_axes y x eye_spacing 0.25 width 600 height 576;

