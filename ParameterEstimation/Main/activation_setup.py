__author__ = 'zwan145'

import os


def activation_create_ipacti(TCa, file_name_r, file_name_w):
    command="awk -v param="+"%12.12f" % ( TCa ) +" -f createIPACTI.awk "+file_name_r+" > "+file_name_w
    os.system(command)