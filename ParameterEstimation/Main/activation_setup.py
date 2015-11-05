__author__ = 'zwan145'

import os


def activation_create_ipacti_exclude_apex(TCa, file_name_r, file_name_w):
    command="awk -v param="+"%12.12f" % ( TCa ) +" -f createIPACTI_exclude_apex.awk "+file_name_r+" > "+file_name_w
    os.system(command)


def activation_get(file_name_r):
    with open(file_name_r, 'r') as f:
        data = f.readlines()
    TCa = float(data[15].split()[-1])
    return TCa