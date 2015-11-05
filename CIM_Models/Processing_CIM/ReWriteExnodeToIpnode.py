"""
This file is designed to rewrite CIM prolate-spheroidal models to CMISS rectangular-cartesian
models 
"""

import os, sys
import scipy
from math import pi
from string import strip

# ############################################################################################
def reWriteExnodeToIpnode(current_study_name):
    ####### Check if the model_ipnode directory exits ######
    if not os.path.exists('../Studies/' + current_study_name + '/model_ipnode'):
        print 'Creating directory model_ipnode for study ', current_study_name
        os.mkdir('../Studies/' + current_study_name + '/model_ipnode')

    os.chdir('../Studies/' + current_study_name)
    ## Define the total number of frames
    total_frame = len(os.listdir('model')) / 2
    print 'The total number of frame for current study is ', total_frame

    ## Start to read CIM files
    for i in range(total_frame):
        ## Define the file name to be read
        if os.path.exists('model/' + current_study_name + '_' + str(i + 1) + '.model.exnode'):
            filename_exnode = ('model/' + current_study_name + '_' + str(i + 1) + '.model.exnode')
        else:
            filename_exnode = ('model/' + current_study_name + '_ZJW_' + str(i + 1) + '.model.exnode')
        ## Read in the CIM model parameters
        [focus, Nodal_Param_CIM] = readExnodeFile(filename_exnode)
        print ' **** Finished Reading nodal parameters for frame .......', i + 1
        print focus
        print filename_exnode
        ############### Reorder the nodes ########################
        total_nodes_CMISS_tmp = 40  ## Node 1 and 18 have 4 versions each
        Nodal_Param_CMISS = scipy.zeros((total_nodes_CMISS_tmp, 7), float)
        ## convert mu,theta to degrees
        Nodal_Param_CIM[:, 5] = Nodal_Param_CIM[:, 5] / pi * 180
        Nodal_Param_CIM[:, 6] = Nodal_Param_CIM[:, 6] / pi * 180

        ## Four versions of Node 1
        Nodal_Param_CMISS[0, :] = Nodal_Param_CIM[36, :]  ## Node 1
        Nodal_Param_CMISS[1, :] = Nodal_Param_CIM[39, :]
        Nodal_Param_CMISS[2, :] = Nodal_Param_CIM[38, :]
        Nodal_Param_CMISS[3, :] = Nodal_Param_CIM[37, :]
        ## Nodes with no versions
        Nodal_Param_CMISS[4, :] = Nodal_Param_CIM[32, :]  ## Node 2
        Nodal_Param_CMISS[5, :] = Nodal_Param_CIM[35, :]  ## Node 3
        Nodal_Param_CMISS[6, :] = Nodal_Param_CIM[34, :]  ## Node 4
        Nodal_Param_CMISS[7, :] = Nodal_Param_CIM[33, :]  ## Node 5
        Nodal_Param_CMISS[8, :] = Nodal_Param_CIM[28, :]  ## Node 6
        Nodal_Param_CMISS[9, :] = Nodal_Param_CIM[31, :]  ## Node 7
        Nodal_Param_CMISS[10, :] = Nodal_Param_CIM[30, :]  ## Node 8
        Nodal_Param_CMISS[11, :] = Nodal_Param_CIM[29, :]  ## Node 9
        Nodal_Param_CMISS[12, :] = Nodal_Param_CIM[24, :]  ## Node 10
        Nodal_Param_CMISS[13, :] = Nodal_Param_CIM[27, :]  ## Node 11
        Nodal_Param_CMISS[14, :] = Nodal_Param_CIM[26, :]  ## Node 12
        Nodal_Param_CMISS[15, :] = Nodal_Param_CIM[25, :]  ## Node 13
        Nodal_Param_CMISS[16, :] = Nodal_Param_CIM[20, :]  ## Node 14
        Nodal_Param_CMISS[17, :] = Nodal_Param_CIM[23, :]  ## Node 15
        Nodal_Param_CMISS[18, :] = Nodal_Param_CIM[22, :]  ## Node 16
        Nodal_Param_CMISS[19, :] = Nodal_Param_CIM[21, :]  ## Node 17
        ## Four versions of Node 18
        Nodal_Param_CMISS[20, :] = Nodal_Param_CIM[16, :]  ## Node 18
        Nodal_Param_CMISS[21, :] = Nodal_Param_CIM[19, :]
        Nodal_Param_CMISS[22, :] = Nodal_Param_CIM[18, :]
        Nodal_Param_CMISS[23, :] = Nodal_Param_CIM[17, :]
        ## Nodes with no versions
        Nodal_Param_CMISS[24, :] = Nodal_Param_CIM[12, :]  ## Node 19
        Nodal_Param_CMISS[25, :] = Nodal_Param_CIM[15, :]  ## Node 20
        Nodal_Param_CMISS[26, :] = Nodal_Param_CIM[14, :]  ## Node 21
        Nodal_Param_CMISS[27, :] = Nodal_Param_CIM[13, :]  ## Node 22
        Nodal_Param_CMISS[28, :] = Nodal_Param_CIM[8, :]  ## Node 23
        Nodal_Param_CMISS[29, :] = Nodal_Param_CIM[11, :]  ## Node 24
        Nodal_Param_CMISS[30, :] = Nodal_Param_CIM[10, :]  ## Node 25
        Nodal_Param_CMISS[31, :] = Nodal_Param_CIM[9, :]  ## Node 26
        Nodal_Param_CMISS[32, :] = Nodal_Param_CIM[4, :]  ## Node 27
        Nodal_Param_CMISS[33, :] = Nodal_Param_CIM[7, :]  ## Node 28
        Nodal_Param_CMISS[34, :] = Nodal_Param_CIM[6, :]  ## Node 29
        Nodal_Param_CMISS[35, :] = Nodal_Param_CIM[5, :]  ## Node 30
        Nodal_Param_CMISS[36, :] = Nodal_Param_CIM[0, :]  ## Node 31
        Nodal_Param_CMISS[37, :] = Nodal_Param_CIM[3, :]  ## Node 32
        Nodal_Param_CMISS[38, :] = Nodal_Param_CIM[2, :]  ## Node 33
        Nodal_Param_CMISS[39, :] = Nodal_Param_CIM[1, :]  ## Node 34
        print ' **** Finished reorder the nodes for frame  .......', i + 1
        ######################### Reorder the nodes finished ##################


        #################### Write CMISS ipnode files #########################
        ## Define the file name to be written
        filename_ipnode = ('model_ipnode/' + current_study_name + '_' + str(i + 1) + '.ipnode')
        writeCMISSIpnode(filename_ipnode, focus, Nodal_Param_CMISS)
        print ' **** Finished writing ipnode for frame  .......', i + 1

        #######################################################################################
        ################### Write CMISS com file to convert models to RC ######
        fileName = 'convert2rc.com'
        writeComfile(fileName, current_study_name, total_frame)
        print ' **** Finished writing the com files .....'

        ################### Write CMGUI com file to visualise models ######
        fileName = 'ViewCine.com'
        writeCMGUIComfile(fileName, current_study_name, total_frame)
        fileName = 'ViewGMModel.com'
        writeCMGUIComfile_CIM(fileName, current_study_name,
                              total_frame)
#############################################################################################



#############################################################################################
def readExnodeFile(fileName):
    """
    This file is designed to read in the exnode file and extract all nodal information
    """

    try:
        file = open(fileName, 'r')
    except IOError:
        print 'ERROR: readExdata_Error: unable to open', fileName
        return


    ## Read in the header (10 lines)
    no_header_line = 6
    for i in range(0, no_header_line):
        header_line = file.readline()
        if (i == 2):  # This line contains the value for focus
            focus_tmp = header_line.split()[6]
            focus = focus_tmp.strip(',')



    ## Define the total number of nodes
    total_nodes_CIM = 40

    ## Define an array
    ##node number,lamda,3 derivatives, mu and theta
    node_param_CIM = scipy.zeros((total_nodes_CIM, 7), float)


    # Start to read the nodal information
    node_info = file.readline()  # Read in the first node

    ## Continue to read until it reaches the end of file
    for i in range(total_nodes_CIM):
        node_info = file.readline()
        node_param_CIM[i, 0] = i + 1
        node_param_CIM[i, 1:5] = ([float(node) for node in node_info.split()[0:4]])  ## lamda and derivatives
        node_info = file.readline()
        node_param_CIM[i, 5] = (float(node_info))  # mu
        node_info = file.readline()
        node_param_CIM[i, 6] = (float(node_info))  # theta

        node_info = file.readline()  # Next node number

    node_param_CIM = scipy.array(node_param_CIM)
    return focus, node_param_CIM


#############################################################################################

#############################################################################################
def writeCMISSIpnode(fileName, focus_current, node_param_CMISS):
    """
    This file is designed to write CMISS ipnode files based on the CIM nodal parameters
    """

    try:
        file = open(fileName, 'w')
    except IOError:
        print 'ERROR: writeCMISSIpnode_Error: unable to open', fileName
        return

    ## Write header information
    file.write(' CMISS Version 2.1  ipnode File Version 2\n')
    file.write(' Heading: Elements created in Perl\n')
    file.write(' \n')
    file.write(' Specify the focus position [1.0]:   {}\n'.format(focus_current))
    file.write(' The number of nodes is [  34]:   34\n')
    file.write(' Number of coordinates [3]: 3\n')
    file.write(' Do you want prompting for different versions of nj=1 [N]? Y\n')
    file.write(' Do you want prompting for different versions of nj=2 [N]? Y\n')
    file.write(' Do you want prompting for different versions of nj=3 [N]? Y\n')
    file.write(' The number of derivatives for coordinate 1 is [0]: 3\n')
    file.write(' The number of derivatives for coordinate 2 is [0]: 3\n')
    file.write(' The number of derivatives for coordinate 3 is [0]: 3\n')
    file.write(' \n')

    total_nodes_CMISS = 34
    ## Write out the nodal parameters
    for k in range(total_nodes_CMISS):
        if (k == 0):
            file.write(' Node number [    1]:     1\n')
            file.write(' The number of versions for nj=1 is [1]: 4\n')
            file.write(' For version number 1:\n')
            file.write(' The Xj(1) coordinate is [ 0.15000E+01]:    {}\n'.format(node_param_CMISS[k, 1]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k, 2]))
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k, 3]))
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k, 4]))
            file.write(' For version number 2:\n')
            file.write(' The Xj(1) coordinate is [ 0.15000E+01]:    {}\n'.format(node_param_CMISS[k + 1, 1]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k + 1, 2]))
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k + 1, 3]))
            file.write(
                ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k + 1, 4]))
            file.write(' For version number 3:\n')
            file.write(' The Xj(1) coordinate is [ 0.15000E+01]:    {}\n'.format(node_param_CMISS[k + 2, 1]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k + 2, 2]))
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k + 2, 3]))
            file.write(
                ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k + 2, 4]))
            file.write(' For version number 4:\n')
            file.write(' The Xj(1) coordinate is [ 0.15000E+01]:    {}\n'.format(node_param_CMISS[k + 3, 1]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k + 3, 2]))
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k + 3, 3]))
            file.write(
                ' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[k + 3, 4]))
            file.write(' The number of versions for nj=2 is [1]: 4\n')
            file.write(' For version number 1:\n')
            file.write(' The Xj(2) coordinate is [ 0.10000E+02]:    {}\n'.format(node_param_CMISS[k, 5]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.000000\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 2:\n')
            file.write(' The Xj(2) coordinate is [ 0.10000E+02]:    {}\n'.format(node_param_CMISS[k, 5]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.000000\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 3:\n')
            file.write(' The Xj(2) coordinate is [ 0.10000E+02]:    {}\n'.format(node_param_CMISS[k, 5]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.000000\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 4:\n')
            file.write(' The Xj(2) coordinate is [ 0.10000E+02]:    {}\n'.format(node_param_CMISS[k, 5]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.000000\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The number of versions for nj=3 is [1]: 4\n')
            file.write(' For version number 1:\n')
            file.write(' The Xj(3) coordinate is [ 0.00000E+00]:     {}\n'.format(node_param_CMISS[k, 6]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 2:\n')
            file.write(' The Xj(3) coordinate is [ 0.00000E+00]:    {}\n'.format(node_param_CMISS[k + 1, 6]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 3:\n')
            file.write(' The Xj(3) coordinate is [ 0.00000E+00]:    {}\n'.format(node_param_CMISS[k + 2, 6]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 4:\n')
            file.write(' The Xj(3) coordinate is [ 0.00000E+00]:    {}\n'.format(node_param_CMISS[k + 3, 6]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' \n')
        elif (k == 17):
            actual_index = k + 3
            file.write(' Node number [    18]:     18\n')
            file.write(' The number of versions for nj=1 is [1]: 4\n')
            file.write(' For version number 1:\n')
            file.write(' The Xj(1) coordinate is [ 0.15000E+01]:    {}\n'.format(node_param_CMISS[actual_index, 1]))
            file.write(
                ' The derivative wrt direction 1 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[actual_index, 2]))
            file.write(
                ' The derivative wrt direction 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[actual_index, 3]))
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: {}\n'.format(
                node_param_CMISS[actual_index, 4]))
            file.write(' For version number 2:\n')
            file.write(' The Xj(1) coordinate is [ 0.15000E+01]:    {}\n'.format(node_param_CMISS[actual_index + 1, 1]))
            file.write(
                ' The derivative wrt direction 1 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[actual_index + 1, 2]))
            file.write(
                ' The derivative wrt direction 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[actual_index + 1, 3]))
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: {}\n'.format(
                node_param_CMISS[actual_index + 1, 4]))
            file.write(' For version number 3:\n')
            file.write(' The Xj(1) coordinate is [ 0.15000E+01]:    {}\n'.format(node_param_CMISS[actual_index + 2, 1]))
            file.write(
                ' The derivative wrt direction 1 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[actual_index + 2, 2]))
            file.write(
                ' The derivative wrt direction 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[actual_index + 2, 3]))
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: {}\n'.format(
                node_param_CMISS[actual_index + 2, 4]))
            file.write(' For version number 4:\n')
            file.write(' The Xj(1) coordinate is [ 0.15000E+01]:    {}\n'.format(node_param_CMISS[actual_index + 3, 1]))
            file.write(
                ' The derivative wrt direction 1 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[actual_index + 3, 2]))
            file.write(
                ' The derivative wrt direction 2 is [ 0.00000E+00]: {}\n'.format(node_param_CMISS[actual_index + 3, 3]))
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: {}\n'.format(
                node_param_CMISS[actual_index + 3, 4]))
            file.write(' The number of versions for nj=2 is [1]: 4\n')
            file.write(' For version number 1:\n')
            file.write(' The Xj(2) coordinate is [ 0.10000E+02]:    {}\n'.format(node_param_CMISS[actual_index, 5]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.000000\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 2:\n')
            file.write(' The Xj(2) coordinate is [ 0.10000E+02]:    {}\n'.format(node_param_CMISS[actual_index, 5]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.000000\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 3:\n')
            file.write(' The Xj(2) coordinate is [ 0.10000E+02]:    {}\n'.format(node_param_CMISS[actual_index, 5]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.000000\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 4:\n')
            file.write(' The Xj(2) coordinate is [ 0.10000E+02]:    {}\n'.format(node_param_CMISS[actual_index, 5]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.000000\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The number of versions for nj=3 is [1]: 4\n')
            file.write(' For version number 1:\n')
            file.write(' The Xj(3) coordinate is [ 0.00000E+00]:     {}\n'.format(node_param_CMISS[actual_index, 6]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 2:\n')
            file.write(' The Xj(3) coordinate is [ 0.00000E+00]:    {}\n'.format(node_param_CMISS[actual_index + 1, 6]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 3:\n')
            file.write(' The Xj(3) coordinate is [ 0.00000E+00]:    {}\n'.format(node_param_CMISS[actual_index + 2, 6]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' For version number 4:\n')
            file.write(' The Xj(3) coordinate is [ 0.00000E+00]:    {}\n'.format(node_param_CMISS[actual_index + 3, 6]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: 0.00000000000000000E+00\n')
            file.write(' \n')
        else:
            file.write(' Node number [    {}]:     {}\n'.format(k + 1, k + 1))
            if (k < 17):
                actual_index = k + 3
            else:
                actual_index = k + 6

            file.write(' The number of versions for nj=1 is [1]: 1\n')
            file.write(' The Xj(1) coordinate is [ 0.24898E+02]:    {}\n'.format(node_param_CMISS[actual_index, 1]))
            file.write(
                ' The derivative wrt direction 1 is [ 0.00000E+00]:    {}\n'.format(node_param_CMISS[actual_index, 2]))
            file.write(
                ' The derivative wrt direction 2 is [ 0.00000E+00]:   {}\n'.format(node_param_CMISS[actual_index, 3]))
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:    {}\n'.format(
                node_param_CMISS[actual_index, 4]))
            file.write(' The number of versions for nj=2 is [1]: 1\n')
            file.write(' The Xj(2) coordinate is [ 0.00000E+00]:    {}\n'.format(node_param_CMISS[actual_index, 5]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]:    0.0000000000000000\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]:    0.0000000000000000\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:    0.0000000000000000\n')
            file.write(' The number of versions for nj=3 is [1]:  1\n')
            file.write(' The Xj(3) coordinate is [ 0.00000E+00]:    {}\n'.format(node_param_CMISS[actual_index, 6]))
            file.write(' The derivative wrt direction 1 is [ 0.00000E+00]:   0.0000000000000000\n')
            file.write(' The derivative wrt direction 2 is [ 0.00000E+00]:    0.0000000000000000\n')
            file.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:   0.0000000000000000\n')
            file.write(' \n')
    file.close()


###############################################################################################################

#######################################################################################
def writeComfile(fileName, current_study_name, total_no_frame):
    """
    This file is written to write the CMISS com file to convert
    """

    try:
        file = open(fileName, 'w')
    except IOError:
        print 'ERROR: writeComfile_Error: unable to open', fileName
        return

    file.write('## This cmiss file converts prolate models to RC models\n')
    file.write('\n')
    file.write('fem define para;r;small;\n')
    file.write('fem define coor;r;prolate;\n')
    file.write('fem define base;r;LVBasis_BiCubicLinear;\n')
    file.write('\n')
    file.write('\n')
    file.write('\n')

    file.write('for ($frame={};$frame<={};$frame++)\n'.format(1, total_no_frame))
    file.write('{\n')
    file.write('\n')
    file.write('	print " ========= Converting Frame $frame ============"; \n')
    file.write('\n')
    file.write('	## Define the file to be read\n')
    str = '	$file_read="model_ipnode/' + current_study_name + '_".$frame;\n'
    file.write(str)
    file.write('	fem define node;r;$file_read;\n')
    file.write('	fem define elem;r;prolate;\n')
    file.write('\n')
    file.write('	## Update Xi1 and Xi2 derivatives for mu (coordinate:2) and theta (coordinate: 3) only\n')
    file.write('	## Note that because there are two nodes with versions, the update needs to take care of the \n')
    file.write('	## versions as well, hence "versions individual"\n')
    file.write('	## mu\n')
    file.write('	fem update nodes derivative 1 linear in 2 wrt xi;\n')
    file.write('	fem update nodes derivative 2 linear in 2 wrt xi versions individual;\n')
    file.write('	## theta\n')
    file.write('	fem update nodes derivative 1 linear in 3 wrt xi;\n')
    file.write('	fem update nodes derivative 2 linear in 3 wrt xi versions individual;\n')
    file.write('\n')
    file.write('	## Write oute the model after updating derivatives\n')
    str = '	$file_read2="model_ipnode/' + current_study_name + '_UpdateDeriv_".$frame;\n'
    file.write(str)
    file.write('	fem define node;w;$file_read2;\n')
    file.write('	## Export the prolate model for visualisation\n')
    file.write('	fem export node;$file_read as prolate;\n')
    file.write('	fem export elem;$file_read as prolate;\n')
    file.write('\n')
    file.write('	## Change coordinate system from 3(prolate) to 1 (RC)\n')
    file.write('	fem change coordinates to 1;\n')
    file.write('\n')
    file.write('	## Write out the converted model\n')
    str = '	$file_write="model_ipnode/' + current_study_name + '_RC_".$frame;\n'
    file.write(str)
    file.write('	fem define nodes;w;$file_write;\n')
    file.write('	fem define element;w;$file_write;\n')
    file.write('	## Export the coverted models for visualisation\n')
    file.write('	fem export node;$file_write as RCAllPhase;\n')
    file.write('	fem export elem;$file_write as RCAllPhase;\n')
    file.write('\n')
    file.write('	## Write in a new ipbase file with tri-cubic interpolation\n')
    file.write('	fem define base;r;LVBasis_Cubic\n')
    file.write('	fem define node;r;$file_write\n')
    file.write('	fem define elem;r;LVModel_Cubic\n')
    file.write('\n')
    file.write('	fem update nodes derivative 3 versions individual;\n')
    file.write('	fem update scale_factor normalise;\n')
    file.write('\n')
    str = '	$file_write="model_ipnode/' + current_study_name + '_RC_Cubic_".$frame;\n'
    file.write(str)
    file.write('	fem define nodes;w;$file_write;\n')
    file.write('	fem define element;w;$file_write;\n')
    file.write('	fem define node;w;temp;\n')
    file.write('	\n');
    file.write('	## Correct xi 2 derivatives\n')
    file.write('	system("perl CorrectXi2Derivatives.sh");\n')
    file.write('	\n')
    file.write('	fem define nodes;r;temp\n')
    file.write('	system("rm temp.ipnode")\n')
    file.write('	fem define nodes;w;$file_write;\n')
    file.write('	fem export node;$file_write as RCAllPhase_Cubic;\n')
    file.write('	fem export elem;$file_write as RCAllPhase_Cubic;\n')
    file.write('\n')
    file.write('	print " ========= Exporting Frame $frame ============";\n')
    file.write('\n')
    file.write('	fem define coor;r;prolate;\n')
    file.write('	fem define base;r;LVBasis_BiCubicLinear;\n')
    file.write('\n')
    file.write('}')
    file.write('\n')
    file.write('fem quit\n')


#######################################################################################

#######################################################################################
def writeCMGUIComfile(fileName, current_study_name, total_no_frame):
    """
    This file is designed to write the cmgui visualisation com file
    """

    try:
        file = open(fileName, 'w')
    except IOError:
        print 'ERROR: writeCMGUIComfile_Error: unable to open', fileName
        return

    file.write('# This CMGUI program is designed to visualize converted CMISS models \n')
    file.write('\n')
    file.write('# Create one window to display the model\n')
    file.write('\n')
    file.write('\n')
    file.write('for ($i={};$i<={};$i++)\n'.format(1, total_no_frame))
    file.write('{\n')
    str = '	$file="model_ipnode/' + current_study_name + '_RC_Cubic_".$i;\n'
    file.write(str)
    file.write('	gfx read node $file time $i;\n')
    file.write('}\n')

    file.write('gfx read element $file;\n')
    file.write(
        'gfx modify g_element RCAllPhase_Cubic general clear circle_discretization 6 default_coordinate coordinates element_discretization "12*12*12" native_discretization none;\n')
    file.write(
        'gfx modify g_element RCAllPhase_Cubic node_points glyph sphere general size "3*3*3" centre 0,0,0 font default select_on material green selected_material default_selected;\n')
    file.write(
        'gfx modify g_element RCAllPhase_Cubic cylinders constant_radius 0.5 select_on material green selected_material default_selected render_shaded;\n')
    file.write(
        'gfx modify g_element RCAllPhase_Cubic surfaces face xi3_0 select_on material green selected_material default_selected render_shaded;\n')

    file.write('# Add the scene editor\n')
    file.write('gfx edit scene;\n')
    file.write('\n')

    file.write('gfx create window 1 double_buffer;\n')
    file.write('gfx modify window 1 image scene default light_model default;\n')
    file.write('gfx modify window 1 image add_light default;\n')
    file.write('gfx modify window 1 layout simple ortho_axes -x -z eye_spacing 0.25 width 1473 height 736;\n')
    file.write('gfx modify window 1 set current_pane 1;\n')
    file.write('gfx modify window 1 background colour 1 1 1 texture none;\n')
    file.write(
        'gfx modify window 1 view parallel eye_point 16.1473 6.22437 -280.88 interest_point 16.1473 6.22437 -1.23459 up_vector -1 0 0 view_angle 45.0019 near_clipping_plane 2.79645 far_clipping_plane 999.357 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;\n')
    file.write('gfx modify window 1 overlay scene none;\n')
    file.write(
        'gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;\n')


#######################################################################################


#######################################################################################
def writeCMGUIComfile_CIM(fileName, current_study_name, total_no_frame):
    """
    This file is designed to write the cmgui visualisation com file
    """

    try:
        file = open(fileName, 'w')
    except IOError:
        print 'ERROR: writeCMGUIComfile_Error: unable to open', fileName
        return

    file.write('# This CMGUI program is designed to visualize converted CMISS models \n')
    file.write('\n')
    file.write('# Create one window to display the model\n')
    file.write('\n')
    file.write('\n')
    file.write('for ($i={};$i<={};$i++)\n'.format(1, total_no_frame))
    file.write('{\n')
    str = '	$file="model/' + current_study_name + '_$i.model.exnode";\n'
    file.write(str)
    file.write('	gfx read node $file time $i;\n')
    file.write('}\n')

    file.write('gfx read element GlobalHermiteParam;\n')
    file.write(
        'gfx modify g_element heart general clear circle_discretization 6 default_coordinate coordinates element_discretization "12*12*12" native_discretization none;\n')
    file.write(
        'gfx modify g_element heart node_points glyph sphere general size "3*3*3" centre 0,0,0 font default select_on material green selected_material default_selected;\n')
    file.write(
        'gfx modify g_element heart cylinders constant_radius 0.5 select_on material green selected_material default_selected render_shaded;\n')
    file.write(
        'gfx modify g_element heart surfaces face xi3_0 select_on material green selected_material default_selected render_shaded;\n')

    file.write('# Add the scene editor\n')
    file.write('gfx edit scene;\n')
    file.write('\n')

    file.write('gfx create window 1 double_buffer;\n')
    file.write('gfx modify window 1 image scene default light_model default;\n')
    file.write('gfx modify window 1 image add_light default;\n')
    file.write('gfx modify window 1 layout simple ortho_axes -x -z eye_spacing 0.25 width 1473 height 736;\n')
    file.write('gfx modify window 1 set current_pane 1;\n')
    file.write('gfx modify window 1 background colour 1 1 1 texture none;\n')
    file.write(
        'gfx modify window 1 view parallel eye_point 16.1473 6.22437 -280.88 interest_point 16.1473 6.22437 -1.23459 up_vector -1 0 0 view_angle 45.0019 near_clipping_plane 2.79645 far_clipping_plane 999.357 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;\n')
    file.write('gfx modify window 1 overlay scene none;\n')
    file.write(
        'gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;\n')


