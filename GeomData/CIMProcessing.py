"""
This code collects the DS model from CIM models and collects the surface points also 
and stores them in the correct folders (diastolic or systolic). 

Author: Jenny Zhinuo Wang
Date: 23rd April 2013

"""
import dircache
import os, sys
import re
import scipy
from shutil import copy

# #######################################################################################

def createSurfacePoints(current_study_name, current_study_frame):
    """
    This file is designed to collect the DS model for all cases and generate
    surface points from all frames and store them in correct Diastolic or
    Systolic folders.
    """

    ############ Step 0: Housekeeping ########################################
    ## Assign variables to the frames
    DS_frame, ED_frame, ES_frame, tot_frame = tuple(current_study_frame)
    DS_frame = int(DS_frame)
    ED_frame = int(ED_frame)
    ES_frame = int(ES_frame)
    tot_frame = int(tot_frame)

    ## Make a directory for the current study
    if not os.path.exists(current_study_name):
        os.mkdir(current_study_name)

    ############# Step 1: Copy the ipnode and ipelem data files ###################
    os.chdir(current_study_name)
    os.system('rm *.*')
    ## Get the current directory
    dst = os.getcwd()
    ## Change to the folder where CMISS ipnode files are stored
    CIMModels_directory = ('../../CIM_Models/Studies/'+current_study_name+'/model_ipnode/')
    os.chdir(CIMModels_directory)
    CIMModels_directory = os.getcwd()

    ## Copy all ipnode and ipelem files into SurfacePoints study folder.
    all_files = os.listdir(CIMModels_directory)
    numIpnode = 1
    numIpelem = 1
    for file in all_files:
        filename, ext = os.path.splitext(file)
        if '_RC_Cubic' in filename:
            if ext == '.ipnode':
                a = re.findall(r'\d+', filename)[1]
                dstfilename = dst + '/' + current_study_name + '_' + str(a) + '.ipnode'
                copy(file, dstfilename)
                numIpnode += 1
            elif ext == '.ipelem':
                a = re.findall(r'\d+', filename)[1]
                dstfilename = dst + '/' + current_study_name + '_' + str(a) + '.ipelem'
                copy(file, dstfilename)
                numIpelem += 1

    ############## Step 2: Sort into Diastole and Systole folders #################

    ################## Active Frames #####################

    os.chdir(dst)
    if not os.path.exists('Active'):
        os.mkdir('Active')
    os.chdir('Active/')
    os.system('rm *.*')
    dstActive = os.getcwd()

    i = ED_frame
    while (i != DS_frame + 1):
        if i > tot_frame:
            i = 1
        filename = dst + '/' + current_study_name + '_' + str(i) + '.ipnode'
        copy(filename, dstActive)
        filename = dst + '/' + current_study_name + '_' + str(i) + '.ipelem'
        copy(filename, dstActive)
        i = i + 1

    ## Auxiliry CMISS files to generate surface points
    all_files = os.listdir('../../TemplateFiles_Surface')
    os.chdir('../../TemplateFiles_Surface')
    for file in all_files:
        copy(file, dstActive)
    os.chdir(dstActive)

    filename = 'CreateSurfacePoints.com'
    i = ED_frame
    while (i != DS_frame + 1):
        if i > tot_frame:
            i = 1
        writeCMISSComFile(filename, current_study_name, str(i))
        print 'Creating surface points for frame ', str(i), ' in active section'
        os.system('cm CreateSurfacePoints.com')
        i = i + 1
    ###############  Passive Frames #############################
    os.chdir(dst)
    if not os.path.exists('Passive'):
        os.mkdir('Passive')
    os.chdir('Passive/')
    os.system('rm *.*')
    dstPassive = os.getcwd()
    print DS_frame
    print tot_frame

    i = DS_frame - 1
    while (i != ED_frame):
        if i + 1 > tot_frame:
            i = 0
        filename = dst + '/' + current_study_name + '_' + str(i + 1) + '.ipnode'
        copy(filename, dstPassive)
        filename = dst + '/' + current_study_name + '_' + str(i + 1) + '.ipelem'
        copy(filename, dstPassive)
        i = i + 1

    ## Auxiliry CMISS files to generate surface points
    all_files = os.listdir('../../TemplateFiles_Surface')
    os.chdir('../../TemplateFiles_Surface')
    for file in all_files:
        copy(file, dstPassive)
    os.chdir(dstPassive)
    filename = 'CreateSurfacePoints.com'
    i = DS_frame - 1
    while (i != ED_frame):
        if i + 1 > tot_frame:
            i = 0
        writeCMISSComFile(filename, current_study_name, str(i + 1))
        print 'Creating surface points for frame ', str(i + 1), ' in passive section'
        os.system('cm CreateSurfacePoints.com')
        i = i + 1

    ################# Step 3: Reset housekeeping ###########################
    os.chdir('../../')


########################################################################################

########################################################################################

def writeCMISSComFile(fileName, current_study_name, current_frame):
    """
    This file is designed to write the CMISS file to generate surface points
    """

    try:
        file = open(fileName, 'w')
    except IOError:
        print 'ERROR: writeCMISSComFile: unable to open', fileName
        return

    file.write(' \n')
    file.write('#### Read in the set-up files\n')
    file.write('fem define para;r;small;			# Define the parameter file\n')
    file.write('fem define coor;r;mapping;			# Define the coordinate file which involves mapping\n')
    file.write('fem define base;r;LVBasisV2\n')
    file.write(' \n')
    file.write('#### Read in the geometric model \n')
    file_str = '$file=' + current_study_name + '_' + current_frame + '\n'
    file.write(file_str)
    file.write('fem define node;r;$file\n')
    file.write('fem define elem;r;$file\n')
    file_str = '$file=' + current_study_name + '_' + current_frame + '\n'
    file.write(file_str)
    file.write('fem define node;w;$file\n')
    file.write('fem define elem;w;$file\n')
    file.write('fem export node;$file as ModelDiastole\n')
    file.write('fem export elem;$file as ModelDiastole\n')
    file.write(' \n')
    file.write('#### Read in the local FE coordinates for endocardial surface #######\n')
    file.write('fem define data;r;Temp			## This is just a temporary data file\n')
    file.write('fem define xi;r;Surface_Points_Endo\n')
    file.write('fem define data;c from_xi\n')
    file_str = '$file=' + current_study_name + '_Surface_Points_Endo_' + current_frame + '\n'
    file.write(file_str)
    file.write('fem define data;w;$file\n')
    file.write('fem list data statistics\n')
    file.write('fem export data;$file as SurfaceEndoDiastole\n')
    file.write('\n')
    file.write('\n')
    file.write('#### Read in the local FE coordinates for epicardial surface #######\n')
    file.write('fem define data;r;Temp			## This is just a temporary data file\n')
    file.write('fem define xi;r;Surface_Points_Epi\n')
    file.write('fem define data;c from_xi\n')
    file_str = '$file=' + current_study_name + '_Surface_Points_Epi_' + current_frame + '\n'
    file.write(file_str)
    file.write('fem define data;w;$file\n')
    file.write('fem list data statistics\n')
    file.write('fem export data;$file as SurfaceEpiDiastole\n')
    file.write('\n')
    file.write('\n')
    file.write('fem reallocate\n')
    file.write('\n')
    file.write('\n')
    file.write('fem quit\n')


########################################################################################
def processPressure(i):
    TOL = 1;
    file = fopen('LV_tracing_' + current_study_name + '.txt')
    pressure_info = file.readline()
    pressure_info = file.readline()
    pressure_info = file.readline()
    i = 0;
    while len(pressure_info) != 0:
        P[i] = pressure_info.split()[0]
        i = i + 1


############################### Main Code #####################################
## Housekeeping
"""
### Set up output log file. 
so = se = open("CreatingSurfaceData.log", 'w', 0)
# re-open stdout without buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

# redirect stdout and stderr to the log file opened above
os.dup2(so.fileno(), sys.stdout.fileno())
os.dup2(se.fileno(), sys.stderr.fileno())
"""
os.system('rm *.*~')

## Extract the important frame numbers for all studies
file = open(os.environ['PARAM_ESTIMATION'] + 'NYStFranFrameNumber.txt', 'r')

study_ID = []
study_frame_tmp = scipy.zeros((0, 3), int)
study_frame = []

study_infor = file.readline()
while len(study_infor) != 0:  # Reaching the end of the file
    a = study_infor.split()[0]
    study_ID.append(a)
    study_frame_tmp = study_infor.split()[1:5]
    study_frame.append(study_frame_tmp)
    no_frames = study_infor.split()[4]
    study_infor = file.readline()

# Calculate the total number of files
no_studies = len(study_ID)

print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print '      The total number of studies is', no_studies
print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

######################### Process all studies ###################################
for i in range(no_studies-3, no_studies-2):
#for i in range(no_studies):
    current_study_name = study_ID[i]
    if current_study_name.find('MR')!= -1:
	    current_study_frame = study_frame[i]
	    print
	    print ''
	    print '*****************************************************************'
	    print '          Current Study Name is ', current_study_name
	    print '*****************************************************************'
	    print ''

	    createSurfacePoints(current_study_name, current_study_frame)

