fem def para;r;small
fem def coor 3,1
fem def base;r;LV_Cubic

fem def node;r;DSWall
fem def elem;r;DSWall
########################################################################################################################
fem def data;r;TEMPLATE_EPI
fem exp data;TEMPLATE_DATA_EPI as Epi_data
fem group faces 5,9,13,16,21,25,29,32,37,41,45,48,53,57,61,64 as EPI
fem def xi;c closest_face faces EPI search_start 20 cross_boundaries
fem list data;TEMPLATE_ERROR_EPI error
fem exp data;TEMPLATE_EX_ERROR_EPI as Epi_error error

########################################################################################################################
fem def data;r;TEMPLATE_ENDO
fem exp data;TEMPLATE_DATA_ENDO as Endo_data
fem group faces 4,8,12,15,20,24,28,31,36,40,44,47,52,56,60,63 as ENDO
fem def xi;c closest_face faces ENDO search_start 20 cross_boundaries
fem list data;TEMPLATE_ERROR_ENDO error
fem exp data;TEMPLATE_EX_ERROR_ENDO as Endo_error error

fem quit
