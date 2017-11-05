__author__ = 'zwan145'

import os


def material_create_ipmate(C1, file_name_r, file_name_w):
    command = "awk -v param=" + "%12.12f" % C1 + " -f createIPMATE.awk "+file_name_r+" > "+file_name_w
    os.system(command)


def material_get(file_name_r):
    with open(file_name_r, 'r') as f:
        data = f.readlines()
    C1 = data[40].split()[-1]
    C1 = float(C1)
    return C1
