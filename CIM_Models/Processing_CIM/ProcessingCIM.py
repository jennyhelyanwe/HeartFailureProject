"""
This python code is designed to process the CIM guide-point models
into a format which can be used in CMISS. 
Author: Vicky Wang

Modified by Jenny Zhinuo Wang 12/02/2014.
"""

import dircache
import os, sys
import shutil
from shutil import copy
import re
from ReWriteExnodeToIpnode import reWriteExnodeToIpnode

# ############################################################################################
def convertCIMModel(current_study_name):
    ## Copy template files
    all_files = os.listdir('TemplateFiles_CIM')
    os.chdir('TemplateFiles_CIM')
    dst = '../../Studies/' + current_study_name
    for file in all_files:
        copy(file, dst)
    ## Change to the Studies directory.
    os.chdir('../../Studies/')

    #################################################################################
    ##### Step 1: Create a folder called model, and copy over all .model files ######
    os.chdir(current_study_name)
    if not os.path.exists('model'):
        os.mkdir('model')

    os.chdir('model_' + current_study_name + '_AY')
    os.system('cp *.model ../model')

    ##### Step 2: Call Jae-Doe's python code to convert .model to exnode files ####
    os.chdir('../')
    os.system('python CIMModelToCAPModel.py')

    print ''
    print '******************************************************************'
    print '  Finished Converting .model to .exnode ', current_study_name
    print '******************************************************************'
    print ''


    ##### Step 3: Convert CIM prolate-spheroidal models to CMISS PS models ########
    os.chdir('../../Processing_CIM/')
    reWriteExnodeToIpnode(current_study_name)

    ##### Step 4: Convert CMISS PS models to CMISS RC models ######################
    os.system('cm convert2rc.com')

    os.chdir('../../Processing_CIM/')

#######################################################################################



############################### Main Python Codes ######################################

os.system('rm *.*~')
## List all files in the CIM_Models folder
all_files_list = dircache.listdir('../Studies/')
print all_files_list
## Work out how many study cases.  
no_studies = len(all_files_list)

print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print '      The total number of studies is', no_studies
print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

for i in range(1,2):
    #for i in range(no_studies):
    current_study_name = all_files_list[i]
	
    print ''
    print '*****************************************************************'
    print '          Current Study Name is ', current_study_name
    print '*****************************************************************'
    print ''

    convertCIMModel(current_study_name)

