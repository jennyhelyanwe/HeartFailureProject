__author__ = 'zwan145'

import os
import re
import numpy
from shutil import copy
from scipy import array
from util_surface_data import readIpdata
import time

# This python script contains functions which deal with pressure and displacement boundary conditions for mechanical
# simulation.


def bc_pressure_get(study_id, study_frame):
    ds, ed, es, tot = tuple(study_frame)
    # This function gets the pressure at current frame from a specified text file.
    f = open(os.environ['HAEMO_DATA'] + study_id + '_registered_LVP_mean.txt', 'r')
    pressure = []
    info = f.readline()
    while len(info) != 0:
        pressure.append(float(info.split()[1]))
        info = f.readline()

    # Offset pressure with DS.
    ds_p = float(pressure[int(ds) - 1])
    for i in range(0, len(pressure)):
        pressure[i] = float(pressure[i]) - ds_p
    return pressure


def bc_pressure_set(p, file_name_r, file_name_w):
    # This function set the pressure in the ipinit file for mechanical simulation.
    try:
        f_r = open(file_name_r, 'r')
        f_w = open(file_name_w, 'w')
    except IOError:
        print 'ERROR: bc_pressure_set: unable to open ', file_name_r
        return

    temp = f_r.readline()
    while temp != '   (3) Restart from previous solution\n':
        temp = f_r.readline()
    temp = f_r.readline()
    toggle = temp.split()[0]
    f_r.close()

    if toggle == '1':
        f_r = open(file_name_r, 'r')
        temp = f_r.readline()
        while temp != ' Do you want to prescribe auxiliary variable/rhs number 1 [N]? Y\n':
            f_w.write(temp)
            temp = f_r.readline()
        f_w.write(temp)
        f_w.write(' The increment is [0.0]:   ' + str(p) + '\n')
        temp = f_r.readline()
        temp = f_r.readline()
        while temp != '':
            f_w.write(temp)
            temp = f_r.readline()
    elif toggle == '2':
        f_r = open(file_name_r, 'r')
        temp = f_r.readline()
        while temp != ' Dependent variable/equation number 4 :\n':
            f_w.write(temp)
            temp = f_r.readline()
        f_w.write(temp)
        temp = f_r.readline()
        f_w.write(temp)
        for i in range(0, 16):
            temp = f_r.readline()
            f_w.write(temp)
            temp = f_r.readline()
            f_w.write(temp)
            junk = f_r.readline()
            f_w.write(' The increment is [0.0]:  ' + str(p) + '\n')
            temp = f_r.readline()
            f_w.write(temp)
        temp = f_r.readline()
        while temp != '':
            f_w.write(temp)
            temp = f_r.readline()
    f_w.close()
    f_r.close()


def bc_displacement_get(study_id):
    # This function gets the closest data points to the four basal nodes for describing the basal displacement
    # boundary condition.
    os.chdir(os.environ['STUDIES'] + study_id + '/LVMechanics' + study_id + '/PassiveMechanics/')

    data_idx = extract_nodal_indices('DSWall.ipnode', 'Surface_Points_Epi_DS.ipdata', 'Surface_Points_Endo_DS.ipdata')
    return data_idx


def bc_displacement_set(node_idx, file_name_current_epi, file_name_next_epi, file_name_current_endo,
                        file_name_next_endo, file_name_rw, file_name_current_model):
    # This function extracts the displacements of the surface data points closest to the basal nodes of the DS model.
    # Extract nodal displacements
    nodes = extract_nodal_disp(node_idx, file_name_current_epi, file_name_next_epi, file_name_current_endo,
                               file_name_next_endo)
    copy(file_name_rw, 'temp.ipinit')

    try:
        f_r = open('temp.ipinit', 'r')
        f_w = open(file_name_rw, 'w')
    except IOError:
        print 'ERROR: bc_displacement_set: unable to open ', file_name_rw
        return

    temp = f_r.readline()
    for component in [0, 1, 2]:
        while temp != ' Enter node #s/name [EXIT]:    31\n':
            f_w.write(temp)
            temp = f_r.readline()
        if component == 0:
            # Match averaged x displacement to avoid apical swing of model.
            print 'LOG: Getting basal node coordinates using '+file_name_current_model
            base_coords = bc_nodal_base_coord_get(file_name_current_model)
            mean_x = extract_ave_x_disp(node_idx, file_name_next_epi, base_coords)
            for i in range(4, 8):
                f_w.write(temp)
                f_w.write(f_r.readline())
                junk = f_r.readline()
                f_w.write(' The increment is [0.0]:    '+str(mean_x)+'\n')
                for j in range(0, 10):
                    f_w.write(f_r.readline())
                temp = f_r.readline()
        else:
            # Match y and z displacements exactly.
            for i in range(4,8):
                f_w.write(temp)
                f_w.write(f_r.readline())
                junk = f_r.readline()
                f_w.write(' The increment is [0.0]:    '+str(nodes[i][component])+'\n')
                print ' The increment is [0.0]:    '+str(nodes[i][component])+'\n'
                for j in range(0, 10):
                    f_w.write(f_r.readline())
                temp = f_r.readline()
    while temp != '':
        f_w.write(temp)
        temp = f_r.readline()

    f_w.close()
    f_r.close()
    print 'LOG: Finished setting basal displacements'


def bc_nodal_base_coord_get(filename):
    ext = filename.split('.')[-1]
    nodes = []
    if ext == 'ipnode':
        with open(filename, 'r') as f:
            temp = f.readline()
            for i in range(0, 1050):
                temp = f.readline()
            for i in range(0, 4):
                coord = []
                for j in range(0, 3):
                    junk = f.readline()
                    temp = f.readline()
                    coord.append(float(temp.split()[-1]))
                    for k in range(0, 7):
                        junk = f.readline()
                junk = f.readline()
                junk = f.readline()
                nodes.append(coord)
            nodes = array(nodes)
    elif ext == 'exnode':
        # Get the coordinates of the basal nodes.
        with open(filename, 'r') as f:
            temp = f.readline()
            fields = f.readline()
            while temp != ' Node:           31\n':
                temp = f.readline()
            for i in range(0, 4):
                coord = []
                for j in range(0, 3):
                    temp = f.readline()
                    coord.append(float(temp.split()[0]))
                    temp = f.readline()
                if fields == ' #Fields=1\n':
                    junk = f.readline()
                else:
                    for i in range(0, 8):
                        junk = f.readline()
                nodes.append(coord)
        nodes = array(nodes)
    return nodes


def extract_ave_x_disp(Node_Data, fileNameNext_Epi, current_base_coords):
    CAPEpi_def = readIpdata(fileNameNext_Epi)
    # Calculate the nodal displacement during inflation
    Node31_x_def = CAPEpi_def[Node_Data[0], 0] - current_base_coords[0, 0] # x displacement measured against model.
    Node32_x_def = CAPEpi_def[Node_Data[1], 0] - current_base_coords[1, 0]
    Node33_x_def = CAPEpi_def[Node_Data[2], 0] - current_base_coords[2, 0]
    Node34_x_def = CAPEpi_def[Node_Data[3], 0] - current_base_coords[3, 0]

    x_def = array([Node31_x_def, Node32_x_def, Node33_x_def, Node34_x_def])
    mean_x = numpy.mean(x_def)
    print 'LOG: Mean x displacement: ', str(mean_x)
    time.sleep(5)
    return mean_x


def extract_nodal_indices(fileName_DS, fileName_surfaceDS_Epi, fileName_surfaceDS_Endo):
    # This function extracts the indices of the data points which are closest to the four basal epicardial nodes of the
    # DS reference model.
    # ####################### Step 1: Read in the ipnode file ########################
    # Read in the fitted ipnode file

    file = open(fileName_DS, 'r')

    # Read a single line at a time
    no_line = 0
    data = file.readline()
    no_line += 1
    Node31 = []
    Node32 = []
    Node33 = []
    Node34 = []
    Node14 = []
    Node15 = []
    Node16 = []
    Node17 = []
    Basal_Node = []

    while no_line <= 1159:
        if no_line == 1053 or no_line == 1062 or no_line == 1071:  # Node 31
            data_infor = re.findall("(-\d+\.\d+|\d+\.\d+|\d+)", data)
            Node31.append(float(data_infor[3]))
            if no_line == 1071:
                #Node31=array(Node31)
                Basal_Node.append(Node31)
        elif no_line == 1082 or no_line == 1091 or no_line == 1100:  # Node 32
            data_infor = re.findall("(-\d+\.\d+|\d+\.\d+|\d+)", data)
            Node32.append(float(data_infor[3]))
            if no_line == 1100:
                #Node32=array(Node32)
                Basal_Node.append(Node32)
        elif no_line == 1111 or no_line == 1120 or no_line == 1129:  # Node 33
            data_infor = re.findall("(-\d+\.\d+|\d+\.\d+|\d+)", data)
            Node33.append(float(data_infor[3]))
            if (no_line == 1129):
                #Node33=array(Node33)
                Basal_Node.append(Node33)
        elif no_line == 1140 or no_line == 1149 or no_line == 1158:  # Node 34
            data_infor = re.findall("(-\d+\.\d+|\d+\.\d+|\d+)", data)
            Node34.append(float(data_infor[3]))
            if no_line == 1158:
                #Node34=array(Node34)
                Basal_Node.append(Node34)
        elif no_line == 476 or no_line == 485 or no_line == 494:      # Node 14
            data_infor = re.findall("(-\d+\.\d+|\d+\.\d+|\d+)", data)
            Node14.append(float(data_infor[3]))
            if no_line == 494:
                #Node31=array(Node31)
                Basal_Node.append(Node14)
        elif no_line == 505 or no_line == 514 or no_line == 523:      # Node 15
            data_infor = re.findall("(-\d+\.\d+|\d+\.\d+|\d+)", data)
            Node15.append(float(data_infor[3]))
            if no_line == 523:
                #Node31=array(Node31)
                Basal_Node.append(Node15)
        elif no_line == 534 or no_line == 543 or no_line == 552:      # Node 16
            data_infor = re.findall("(-\d+\.\d+|\d+\.\d+|\d+)", data)
            Node16.append(float(data_infor[3]))
            if no_line == 552:
                #Node31=array(Node31)
                Basal_Node.append(Node16)
        elif no_line == 563 or no_line == 572 or no_line == 581:      # Node 17
            data_infor = re.findall("(-\d+\.\d+|\d+\.\d+|\d+)", data)
            Node17.append(float(data_infor[3]))
            if no_line == 581:
                #Node31=array(Node31)
                Basal_Node.append(Node17)
        data = file.readline()
        no_line += 1

    # Step 2: read in the data at DS
    Node_Data = []

    CAPEndo = readIpdata(fileName_surfaceDS_Endo)
    no_data = len(CAPEndo)
    # Loop through all data
    for i in range(0,4):
        Basal_Node_tmp = array(Basal_Node[i])
        print Basal_Node_tmp
        Eucli_Disp = []
        for j in range(no_data):
            data = CAPEndo[j]
            # Calculate the euclidean distance
            Eucli_Disp.append(numpy.linalg.norm(data - Basal_Node_tmp))

        # Calculate the minimum
        Node_ED_min = min(Eucli_Disp)
        Node_Data.append(Eucli_Disp.index(Node_ED_min))

    CAPEpi = readIpdata(fileName_surfaceDS_Epi)
    no_data = len(CAPEpi)
    # Loop through all data
    for i in range(4,8):
        Basal_Node_tmp = array(Basal_Node[i])
        Eucli_Disp = []
        for j in range(no_data):
            data = CAPEpi[j]
            # Calculate the euclidean distance
            Eucli_Disp.append(numpy.linalg.norm(data - Basal_Node_tmp))

        # Calculate the minimum
        Node_ED_min = min(Eucli_Disp)
        Node_Data.append(Eucli_Disp.index(Node_ED_min))

    # Print the four data index
    print 'The data indices for Node 14,15,16,17,31,32,32,34 are ', Node_Data
    return Node_Data


def extract_nodal_disp(Node_Data, fileNameCurrent_Epi, fileNameNext_Epi, fileNameCurrent_Endo, fileNameNext_Endo):
    # This function extracts the displacement at the current frame of the data points which are closest to basal
    # epicardial nodes of the DS model at the current frame in the cardiac cycle.
    # Read in the data at ED and ES
    CAPEpi = readIpdata(fileNameCurrent_Epi)
    CAPEpi_def = readIpdata(fileNameNext_Epi)

    print '================ Write out ipinit file based on displacement  ===================='
    # Calculate the nodal displacement during inflation
    Node31_x_def = CAPEpi_def[Node_Data[0], 0] - CAPEpi[Node_Data[0], 0]
    Node32_x_def = CAPEpi_def[Node_Data[1], 0] - CAPEpi[Node_Data[1], 0]
    Node33_x_def = CAPEpi_def[Node_Data[2], 0] - CAPEpi[Node_Data[2], 0]
    Node34_x_def = CAPEpi_def[Node_Data[3], 0] - CAPEpi[Node_Data[3], 0]

    Node31_y_def = CAPEpi_def[Node_Data[0], 1] - CAPEpi[Node_Data[0], 1]
    Node32_y_def = CAPEpi_def[Node_Data[1], 1] - CAPEpi[Node_Data[1], 1]
    Node33_y_def = CAPEpi_def[Node_Data[2], 1] - CAPEpi[Node_Data[2], 1]
    Node34_y_def = CAPEpi_def[Node_Data[3], 1] - CAPEpi[Node_Data[3], 1]

    Node31_z_def = CAPEpi_def[Node_Data[0], 2] - CAPEpi[Node_Data[0], 2]
    Node32_z_def = CAPEpi_def[Node_Data[1], 2] - CAPEpi[Node_Data[1], 2]
    Node33_z_def = CAPEpi_def[Node_Data[2], 2] - CAPEpi[Node_Data[2], 2]
    Node34_z_def = CAPEpi_def[Node_Data[3], 2] - CAPEpi[Node_Data[3], 2]

    print 'Nodal displacement: '
    print 'Node 31 x, y, z diplacements = ', Node31_x_def, Node31_y_def, Node31_z_def
    print 'Node 32 x, y, z diplacements = ', Node32_x_def, Node32_y_def, Node32_z_def
    print 'Node 33 x, y, z diplacements = ', Node33_x_def, Node33_y_def, Node33_z_def
    print 'Node 34 x, y, z diplacements = ', Node34_x_def, Node34_y_def, Node34_z_def

    Node31_def = [Node31_x_def, Node31_y_def, Node31_z_def]
    Node32_def = [Node32_x_def, Node32_y_def, Node32_z_def]
    Node33_def = [Node33_x_def, Node33_y_def, Node33_z_def]
    Node34_def = [Node34_x_def, Node34_y_def, Node34_z_def]

    # Read in the data at ED and ES
    CAPEndo = readIpdata(fileNameCurrent_Endo)
    CAPEndo_def = readIpdata(fileNameNext_Endo)

    print '================ Write out ipinit file based on displacement  ===================='
    # Calculate the nodal displacement during inflation
    Node14_x_def = CAPEndo_def[Node_Data[0], 0] - CAPEndo[Node_Data[0], 0]
    Node15_x_def = CAPEndo_def[Node_Data[1], 0] - CAPEndo[Node_Data[1], 0]
    Node16_x_def = CAPEndo_def[Node_Data[2], 0] - CAPEndo[Node_Data[2], 0]
    Node17_x_def = CAPEndo_def[Node_Data[3], 0] - CAPEndo[Node_Data[3], 0]

    Node14_y_def = CAPEndo_def[Node_Data[0], 1] - CAPEndo[Node_Data[0], 1]
    Node15_y_def = CAPEndo_def[Node_Data[1], 1] - CAPEndo[Node_Data[1], 1]
    Node16_y_def = CAPEndo_def[Node_Data[2], 1] - CAPEndo[Node_Data[2], 1]
    Node17_y_def = CAPEndo_def[Node_Data[3], 1] - CAPEndo[Node_Data[3], 1]

    Node14_z_def = CAPEndo_def[Node_Data[0], 2] - CAPEndo[Node_Data[0], 2]
    Node15_z_def = CAPEndo_def[Node_Data[1], 2] - CAPEndo[Node_Data[1], 2]
    Node16_z_def = CAPEndo_def[Node_Data[2], 2] - CAPEndo[Node_Data[2], 2]
    Node17_z_def = CAPEndo_def[Node_Data[3], 2] - CAPEndo[Node_Data[3], 2]

    print 'Nodal displacement: '
    print 'Node 14 x, y, z diplacements = ', Node14_x_def, Node14_y_def, Node14_z_def
    print 'Node 15 x, y, z diplacements = ', Node15_x_def, Node15_y_def, Node15_z_def
    print 'Node 16 x, y, z diplacements = ', Node16_x_def, Node16_y_def, Node16_z_def
    print 'Node 17 x, y, z diplacements = ', Node17_x_def, Node17_y_def, Node17_z_def

    Node14_def = [Node14_x_def, Node14_y_def, Node14_z_def]
    Node15_def = [Node15_x_def, Node15_y_def, Node15_z_def]
    Node16_def = [Node16_x_def, Node16_y_def, Node16_z_def]
    Node17_def = [Node17_x_def, Node17_y_def, Node17_z_def]

    return Node14_def, Node15_def, Node16_def, Node17_def, Node31_def, Node32_def, Node33_def, Node34_def
