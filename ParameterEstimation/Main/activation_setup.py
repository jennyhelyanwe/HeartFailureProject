__author__ = 'zwan145'

import os


def activation_create_ipacti(TCa, TCa_apex, file_name_r, file_name_w):
    command = "awk -v param=" + "%12.12f" % TCa + " '{if((NR==18)||(NR==22)){$NF=param;printf(\" %s\\n\",$0);}" \
                                                  "else{print}}' "+file_name_r+" > temp.ipacti"
    os.system(command)
    command = "awk -v param="+"%12.12f" % TCa_apex + " '{if(NR==26){$NF=param;printf(\" %s\\n\",$0);}" \
                                                     "else{print}}' temp.ipacti > "+file_name_w
    os.system(command)


def activation_get(file_name_r):
    with open(file_name_r, 'r') as f:
        data = f.readlines()
    TCa = float(data[17].split()[-1])
    return TCa
