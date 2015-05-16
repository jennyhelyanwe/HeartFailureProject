"""
This python code is designed to extract the DS ipnode file for each case and 
generate surface points at ED and ES frames
Author: Vicky Wang
"""

import dircache
import os,sys
import re
import scipy
from shutil import copy

########################################################################################

def creatsSurfacePoints(current_study_name,current_study_frame) :

	"""
	This file is designed to collect the DS model for all cases and generate
	surface points from ED and ES frame
	"""

	## Change directory
	os.chdir('../GenerateSurfacePoints')
	## Make a directory for the current study
	if not os.path.exists (current_study_name) :
		os.mkdir(current_study_name)


	############# Step 1: Copy the ipnode files  #####################
	os.chdir(current_study_name)
	## Get the current directory
	dst=os.getcwd()
	## Change to the folder where CMISS ipnode files are stored
	CMISS_directory=('../../CIMModels/' + current_study_name + '/model_ipnode/')
	os.chdir(CMISS_directory)

	## Assign variables to the frames
	DS_frame,ED_frame,ES_frame=tuple(current_study_frame)
	## Copy the ipnode file at DS, ED and ES frame as well as the ipelem file
	CMISS_Ipnode_DS=current_study_name + '_RC_Cubic_' + str(DS_frame) + '.ipnode'
	copy(CMISS_Ipnode_DS,dst)
	CMISS_Ipnode_ED=current_study_name + '_RC_Cubic_' + str(ED_frame) + '.ipnode'
	copy(CMISS_Ipnode_ED,dst)
	CMISS_Ipnode_ES=current_study_name + '_RC_Cubic_' + str(ES_frame) + '.ipnode'
	copy(CMISS_Ipnode_ES,dst)
	CMISS_Ipelem_DS=current_study_name + '_RC_Cubic_' + str(DS_frame) + '.ipelem'
	copy(CMISS_Ipelem_DS,dst)	
	CMISS_Ipelem_ED=current_study_name + '_RC_Cubic_' + str(ED_frame) + '.ipelem'
	copy(CMISS_Ipelem_ED,dst)	
	CMISS_Ipelem_ES=current_study_name + '_RC_Cubic_' + str(ES_frame) + '.ipelem'
	copy(CMISS_Ipelem_ES,dst)	

	## Change back directory
	os.chdir(dst)
	## Auxiliry CMISS files to generate surface points
	all_files=os.listdir('../TemplateFiles')
	os.chdir('../TemplateFiles')
	for file in all_files :
		copy(file,dst)
	os.chdir(dst)
	
	############ Step 2: Generate surface points ######################
	filename='CreateSurfacePoints.com'
	writeCMISSComFile(filename,current_study_name,DS_frame,ED_frame,ES_frame)
	os.system('cm CreateSurfacePoints.com')

	os.chdir('../../PYTHON_Codes_ProcessingCIM')

########################################################################################

########################################################################################

def writeCMISSComFile(fileName,current_study_name,DS_frame,ED_frame,ES_frame) :

	"""
	This file is designed to write the CMISS file to generate surface points
	"""

	try:
		file = open( fileName, 'w' )
	except IOError:
		print 'ERROR: writeCMISSComFile: unable to open', fileName
		return

	file.write(' \n')	
	file.write('##################################### DS #####################################\n')
	file.write('#### Read in the set-up files\n')
	file.write('fem define para;r;small;			# Define the parameter file\n')
	file.write('fem define coor;r;mapping;			# Define the coordinate file which involves mapping\n')
	file.write('fem define base;r;LVBasisV2\n')
	file.write(' \n')
	file.write('#### Read in the geometric model at DS \n')
	file_str='$file=' + current_study_name + '_RC_Cubic_' + str(DS_frame) + '\n'
	file.write(file_str)
	file.write('fem define node;r;$file\n')
	file.write('fem define elem;r;$file\n')
	file_str='$file=' + current_study_name + '_DS_Model\n'
	file.write(file_str)
	file.write('fem define node;w;$file\n')
	file.write('fem define elem;w;$file\n')
	file.write('fem export node;$file as $file\n')
	file.write('fem export elem;$file as $file\n')
	file.write(' \n')
	file.write('#### Read in the local FE coordinates for endocardial surface #######\n')
	file.write('fem define data;r;Temp			## This is just a temporary data file\n')
	file.write('fem define xi;r;Surface_Points_Endo\n')
	file.write('fem define data;c from_xi\n')
	file_str='$file=' + current_study_name + '_Surface_Points_Endo_DS\n'
	file.write(file_str)
	file.write('fem define data;w;$file\n')
	file.write('fem list data statistics\n')
	file.write('fem export data;$file as $file\n')
	file.write('\n')
	file.write('\n')
	file.write('#### Read in the local FE coordinates for epicardial surface #######\n')
	file.write('fem define data;r;Temp			## This is just a temporary data file\n')
	file.write('fem define xi;r;Surface_Points_Epi\n')
	file.write('fem define data;c from_xi\n')
	file_str='$file=' + current_study_name + '_Surface_Points_Epi_DS\n'
	file.write(file_str)
	file.write('fem define data;w;$file\n')
	file.write('fem list data statistics\n')
	file.write('fem export data;$file as $file\n')
	file.write('\n')
	file.write('\n')
	file.write('fem reallocate\n')
	file.write('\n')
	file.write('##################################### ED #####################################\n')
	file.write('#### Read in the set-up files\n')
	file.write('fem define para;r;small;			# Define the parameter file\n')
	file.write('fem define coor;r;mapping;			# Define the coordinate file which involves mapping\n')
	file.write('fem define base;r;LVBasisV2\n')
	file.write(' \n')
	file.write('#### Read in the geometric model at ED \n')
	file_str='$file=' + current_study_name + '_RC_Cubic_' + str(ED_frame) + '\n'
	file.write(file_str)
	file.write('fem define node;r;$file\n')
	file.write('fem define elem;r;$file\n')
	file_str='$file=' + current_study_name + '_ED_Model\n'
	file.write(file_str)
	file.write('fem define node;w;$file\n')
	file.write('fem define elem;w;$file\n')
	file.write('fem export node;$file as $file\n')
	file.write('fem export elem;$file as $file\n')
	file.write('\n')
	file.write('#### Read in the local FE coordinates for endocardial surface #######\n')
	file.write('fem define data;r;Temp			## This is just a temporary data file\n')
	file.write('fem define xi;r;Surface_Points_Endo\n')
	file.write('fem define data;c from_xi\n')
	file_str='$file=' + current_study_name + '_Surface_Points_Endo_ED\n'
	file.write(file_str)
	file.write('fem define data;w;$file\n')
	file.write('fem list data statistics\n')
	file.write('fem export data;$file as $file\n')
	file.write('\n')
	file.write('\n')
	file.write('#### Read in the local FE coordinates for epicardial surface #######\n')
	file.write('fem define data;r;Temp			## This is just a temporary data file\n')
	file.write('fem define xi;r;Surface_Points_Epi\n')
	file.write('fem define data;c from_xi\n')
	file_str='$file=' + current_study_name + '_Surface_Points_Epi_ED\n'
	file.write(file_str)
	file.write('fem define data;w;$file\n')
	file.write('fem list data statistics\n')
	file.write('fem export data;$file as $file\n')
	file.write('\n')
	file.write('\n')
	file.write('fem reallocate\n')
	file.write('\n')
	file.write('##################################### ES #####################################\n')
	file.write('#### Read in the set-up files\n')
	file.write('fem define para;r;small;			# Define the parameter file\n')
	file.write('fem define coor;r;mapping;			# Define the coordinate file which involves mapping\n')
	file.write('fem define base;r;LVBasisV2\n')
	file.write('\n')
	file.write('#### Read in the geometric model at ES \n')
	file_str='$file=' + current_study_name + '_RC_Cubic_' + str(ES_frame) + '\n'
	file.write(file_str)
	file.write('fem define node;r;$file\n')
	file.write('fem define elem;r;$file\n')
	file_str='$file=' + current_study_name + '_ES_Model\n'
	file.write(file_str)
	file.write('fem define node;w;$file\n')
	file.write('fem define elem;w;$file\n')
	file.write('fem export node;$file as $file\n')
	file.write('fem export elem;$file as $file\n')
	file.write('\n')
	file.write('\n')
	file.write('#### Read in the local FE coordinates for endocardial surface #######\n')
	file.write('fem define data;r;Temp			## This is just a temporary data file\n')
	file.write('fem define xi;r;Surface_Points_Endo\n')
	file.write('fem define data;c from_xi\n')
	file_str='$file=' + current_study_name + '_Surface_Points_Endo_ES\n'
	file.write(file_str)
	file.write('fem define data;w;$file\n')
	file.write('fem list data statistics\n')
	file.write('fem export data;$file as $file\n')
	file.write('\n')
	file.write('\n')
	file.write('#### Read in the local FE coordinates for epicardial surface #######\n')
	file.write('fem define data;r;Temp			## This is just a temporary data file\n')
	file.write('fem define xi;r;Surface_Points_Epi\n')
	file.write('fem define data;c from_xi\n')
	file_str='$file=' + current_study_name + '_Surface_Points_Epi_ES\n'
	file.write(file_str)
	file.write('fem define data;w;$file\n')
	file.write('fem list data statistics\n')
	file.write('fem export data;$file as $file\n')
	file.write('\n')
	file.write('\n')	
	file.write('fem quit\n')

########################################################################################

###############################Main Python Codes ######################################

os.system('rm *.*~')

## Read a txt file which lists the frame number for all studies
file = open( '../GenerateSurfacePoints/NYStFran_FrameNumber.txt', 'r' )

study_ID=[]
study_frame_tmp=scipy.zeros((0,3), int)
study_frame=[]

study_infor=file.readline()

while len(study_infor) != 0 :	# Reaching the end of the file
	
	a=study_infor.split()[0]
	study_ID.append(a)
	study_frame_tmp=study_infor.split()[1:4]
	study_frame.append(study_frame_tmp)
	study_infor=file.readline()

## Calculate the total number of files
no_studies=len(study_ID)


print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print '      The total number of studies is',no_studies 
print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

######################### Statr to process all models ###################################
for i in range(no_studies):
	current_study_name=study_ID[i]
	current_study_frame=study_frame[i]
	
	print 
	print ''
	print '*****************************************************************'
	print '          Current Study Name is ',current_study_name
	print '*****************************************************************'
	print ''
	creatsSurfacePoints(current_study_name,current_study_frame)
	

######################## Categories model and data files ###############################
os.chdir('../GenerateSurfacePoints')
for i in range(no_studies):
	current_study_name=study_ID[i]
	print 
	print ''
	print '*****************************************************************'
	print '          Current Study Name is ',current_study_name
	print '*****************************************************************'
	print ''
	os.chdir(current_study_name)
	os.system("cp `find . -name '*DS_Model.ipnode*'`  ../../Trail/Data")
	os.system("cp `find . -name '*DS_Model.ipelem*'`  ../../Trail/Data")
	os.system("cp `find . -name '*_DS.ipdata*'`  ../../Trail/Data")	
	os.system("cp `find . -name '*_ED.ipdata*'`  ../../Trail/Data")	
	os.system("cp `find . -name '*_ES.ipdata*'`  ../../Trail/Data")
	os.chdir('../')

	os.chdir('../')		

	
	
