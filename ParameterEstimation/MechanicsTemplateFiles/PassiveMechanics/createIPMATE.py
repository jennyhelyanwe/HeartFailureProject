"""
This python script is called from the CMISS command line to rewrite the IPMATE file to the value given in the commandline inputs. 

Author: Zhinuo Jenny Wang
Date: 14th July 2014
"""

import sys
import os


#==============================================================================================================================#
def createIPMATE(filename,passive_parameters,filenamew):
	comand="awk -v param="+"%12.12f" % ( C1 ) +" -f createIPMATE.awk "+filename+" > "+filenamew

	os.system(comand)

input = sys.argv
print input

index = input[1]
print index
filename = 'LV_CubicGuc_TEMPLATE.ipmate'
filenamew = 'output_debug/LV_CubicGuc_'+str(index)+'.ipmate'
print filename
print filenamew

C1 = float(input[2])
print C1

createIPMATE(filename, C1, filenamew)

