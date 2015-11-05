__author__ = 'zwan145'

import os
from shutil import copy
import scipy
from numpy import array

# This python script contains functions for setting up the geometry of the LV wall and cavity models.


def geom_setup_ref_cavity(study_id, study_frame):
    # This function copies over DS reference model into study folder, and generates the LV cavity model.
    # Get frame numbers
    ds, ed, es, tot = tuple(study_frame)

    # Copy over wall geometry into working folder
    dir_from = os.environ['GEOM_DATA'] + study_id + '/Passive/'
    dir_work = os.environ['STUDIES'] + study_id + '/GeometricModel' + study_id

    copy(dir_from + study_id + '_' + str(ds) + '.ipnode', dir_work + '/DSWall.ipnode')
    copy(dir_from + study_id + '_' + str(ds) + '.ipelem', dir_work + '/DSWall.ipelem')

    # Copy template files into working folder
    all_f = os.listdir(os.environ['GEOM_TEMPLATE'])
    os.chdir(os.environ['GEOM_TEMPLATE'])
    for f in all_f:
        copy(f, dir_work)
    os.chdir(dir_work)

    # Create cavity model
    create_cavity_model('DSWall.ipnode', 'DSCavity.ipnode')

    # Export wall and cavity models for visualisation
    os.system('cm LVHumanDSModel.com')
#
#=======================================================================================================================
#


def geom_setup_data(study_id, study_frame):
    # This function copies reference model and important surface data points to mechanics working folder.
    # Get frame numbers
    ds, ed, es, tot = tuple(study_frame)

    # Copy reference model into mechanics folders
    dir_from = os.environ['STUDIES'] + study_id + '/GeometricModel' + study_id
    dir_to_p = os.environ['STUDIES'] + study_id + '/LVMechanics' + study_id + '/PassiveMechanics'
    dir_to_a = os.environ['STUDIES'] + study_id + '/LVMechanics' + study_id + '/ActiveMechanics'

    copy(dir_from+'/DSWall.ipnode', dir_to_p)
    copy(dir_from+'/DSWall.ipelem', dir_to_p)
    copy(dir_from+'/DSWall.ipnode', dir_to_a)
    copy(dir_from+'/DSWall.ipelem', dir_to_a)
    copy(dir_from+'/DSCavity.ipnode', dir_to_p)
    copy(dir_from+'/DSCavity.ipelem', dir_to_p)
    copy(dir_from+'/DSCavity.ipnode', dir_to_a)
    copy(dir_from+'/DSCavity.ipelem', dir_to_a)

    # Copy important surface data frames to mechanics folders.
    dir_from = os.environ['GEOM_DATA']+study_id+'/Passive/'
    copy(dir_from+study_id+'_Surface_Points_Endo_'+str(ds)+'.ipdata', dir_to_p+'/Surface_Points_Endo_DS.ipdata')
    copy(dir_from+study_id+'_Surface_Points_Epi_'+str(ds)+'.ipdata', dir_to_p+'/Surface_Points_Epi_DS.ipdata')
    copy(dir_from+study_id+'_Surface_Points_Endo_'+str(ed)+'.ipdata', dir_to_p+'/Surface_Points_Endo_ED.ipdata')
    copy(dir_from+study_id+'_Surface_Points_Epi_'+str(ed)+'.ipdata', dir_to_p+'/Surface_Points_Epi_ED.ipdata')

    dir_from = os.environ['GEOM_DATA']+study_id+'/Active/'
    copy(dir_from+study_id+'_Surface_Points_Endo_'+str(es)+'.ipdata', dir_to_a+'/Surface_Points_Endo_ES.ipdata')
    copy(dir_from+study_id+'_Surface_Points_Epi_'+str(es)+'.ipdata', dir_to_a+'/Surface_Points_Epi_ES.ipdata')
#
#=======================================================================================================================
#


def create_cavity_model(file_name_wall, file_name_cavity):
    # This function creates the LV cavity model using the endocardial nodes of the LV wall model.
    try:
        file_wall = open(file_name_wall, 'r')
        file_cavity = open(file_name_cavity, 'w')
    except IOError:
        print 'ERROR: CreateCavityModel: unable to open', file_name_wall
        return

    no_line = 0
    data = file_wall.readline()
    no_line = no_line + 1

    ## Initialise a variable to store the x coordinates of the basal nodes
    x_coor_base = []
    x_coor_apex = []
    y_coor_base = []
    z_coor_base = []

    while (no_line <= 588):  ## Line 588 is the end of node 17
        if (no_line == 16):  ## x-coordinate of Node 1
            data_coor_tmp = data.split()
            data_coor = data_coor_tmp[len(data_coor_tmp) - 1]
            x_coor_apex.append(float(data_coor))
        elif (no_line == 476):  ## x-coordinate of Node 14
            data_coor_tmp = data.split()
            data_coor = data_coor_tmp[len(data_coor_tmp) - 1]
            x_coor_base.append(float(data_coor))
        elif (no_line == 485):  ## y-coordinate of Node 14
            data_coor_tmp = data.split()
            data_coor = data_coor_tmp[len(data_coor_tmp) - 1]
            y_coor_base.append(float(data_coor))
        elif (no_line == 505):  ## x-coordinate of Node 15
            data_coor_tmp = data.split()
            data_coor = data_coor_tmp[len(data_coor_tmp) - 1]
            x_coor_base.append(float(data_coor))
        elif (no_line == 523):  ## z-coordinate of Node 15
            data_coor_tmp = data.split()
            data_coor = data_coor_tmp[len(data_coor_tmp) - 1]
            z_coor_base.append(float(data_coor))
        elif (no_line == 534):  ## x-coordinates of Node 16
            data_coor_tmp = data.split()
            data_coor = data_coor_tmp[len(data_coor_tmp) - 1]
            x_coor_base.append(float(data_coor))
        elif (no_line == 543):  ## y-coordinate of Node 16
            data_coor_tmp = data.split()
            data_coor = data_coor_tmp[len(data_coor_tmp) - 1]
            y_coor_base.append(float(data_coor))
        elif (no_line == 563):  ## x-coordinate of Node 17
            data_coor_tmp = data.split()
            data_coor = data_coor_tmp[len(data_coor_tmp) - 1]
            x_coor_base.append(float(data_coor))
        elif (no_line == 581):  ## z-coordinate of Node 17
            data_coor_tmp = data.split()
            data_coor = data_coor_tmp[len(data_coor_tmp) - 1]
            z_coor_base.append(float(data_coor))

        if (no_line == 4):
            string = ' The number of nodes is [    21]:     21\n'
        else:
            string = data

        file_cavity.write(str(string))
        data = file_wall.readline()
        no_line = no_line + 1

    file_wall.close()
    x_coor_apex = array(x_coor_apex)
    x_coor_base = array(x_coor_base)
    y_coor_base = array(y_coor_base)
    z_coor_base = array(z_coor_base)
    print 'The x-coordinate of the apex is ', float(x_coor_apex)
    print 'The x-coordinates of the base is ', x_coor_base
    print 'The y-coordinates of the base is ', y_coor_base
    print 'The z-coordinates of the base is ', z_coor_base
    ## Calculate the mean coordinate of the basal node
    x_coor_mean = scipy.mean(x_coor_base)
    y_coor_mean = scipy.mean(y_coor_base)
    z_coor_mean = scipy.mean(z_coor_base)
    print 'The mean x-coordinate of the base is ', float(x_coor_mean)
    print 'The mean y-coordinate of the base is ', float(y_coor_mean)
    print 'The mean z-coordinate of the base is ', float(z_coor_mean)

    ## Calculate the base-to-apex dimension
    base_apex_dimen = x_coor_apex - x_coor_mean
    print 'The Base-To-Apex dimension is ', float(base_apex_dimen)
    ### Calculate the coordinates of other nodes in the cavity
    x_coor_118 = float(x_coor_apex - 4 * base_apex_dimen / 4)
    x_coor_119 = float(x_coor_apex - 3 * base_apex_dimen / 4)
    x_coor_120 = float(x_coor_apex - 2 * base_apex_dimen / 4)
    x_coor_121 = float(x_coor_apex - 1 * base_apex_dimen / 4)
    print 'The x-coordinate of node 118 is ', x_coor_118
    print 'The x-coordinate of node 119 is ', x_coor_119
    print 'The x-coordinate of node 120 is ', x_coor_120
    print 'The x-coordinate of node 121 is ', x_coor_121

    ## Add these information to the cavity ipnode file
    string = ' \n'
    file_cavity.write(str(string))
    string = ' Node number [  118]:   118\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=1 is [1]:  4\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_118) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_118) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_118) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_118) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=2 is [1]:  4\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                 ' + str(y_coor_mean) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                 ' + str(y_coor_mean) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                 ' + str(y_coor_mean) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                 ' + str(y_coor_mean) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=3 is [1]:  4\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                 ' + str(z_coor_mean) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                 ' + str(z_coor_mean) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                 ' + str(z_coor_mean) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                 ' + str(z_coor_mean) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = '  \n'
    file_cavity.write(str(string))
    string = ' Node number [  119]:   119\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=1 is [1]:  4\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_119) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_119) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_119) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_119) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=2 is [1]:  4\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=3 is [1]:  4\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = '  \n'
    file_cavity.write(str(string))
    string = ' Node number [  120]:   120\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=1 is [1]:  4\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_120) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_120) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_120) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_120) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=2 is [1]:  4\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=3 is [1]:  4\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = '  \n'
    file_cavity.write(str(string))
    string = ' Node number [  121]:   121\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=1 is [1]:  8\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_121) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_121) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_121) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_121) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 5:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_121) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 6:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_121) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 7:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_121) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 8:\n'
    file_cavity.write(str(string))
    string = ' The Xj(1) coordinate is [ 0.00000E+00]:                 ' + str(x_coor_121) + '\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=2 is [1]:  8\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 5:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 6:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 7:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 8:\n'
    file_cavity.write(str(string))
    string = ' The Xj(2) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The number of versions for nj=3 is [1]:  8\n'
    file_cavity.write(str(string))
    string = ' For version number 1:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 2:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 3:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 4:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 5:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 6:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 7:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' For version number 8:\n'
    file_cavity.write(str(string))
    string = ' The Xj(3) coordinate is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 1 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt direction 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))
    string = ' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:                         0\n'
    file_cavity.write(str(string))

    file_cavity.close()

    return
#
#=======================================================================================================================
#