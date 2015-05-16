__author__ = 'zwan145'

import os


def material_create_ipmate(C1, file_name_r, file_name_w):
    command = "awk -v param=" + "%12.12f" % (C1 ) + " -f createIPMATE.awk "+file_name_r+" > "+file_name_w
    os.system(command)