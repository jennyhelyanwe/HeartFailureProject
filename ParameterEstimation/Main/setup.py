__author__ = 'zwan145'

import os
import sys

# This script contains functions used for general program set up.


def setup_log_file(study_id, log_name):
    # This function writes all outputs to log file.
    os.chdir(os.environ['STUDIES']+study_id)
    so = se = open(log_name, 'w', 0)
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    os.dup2(so.fileno(), sys.stdout.fileno())
    os.dup2(se.fileno(), sys.stderr.fileno())
#
#=======================================================================================================================
#


def setup_dir(study_id):
    # This function set up the study directories.
    study_dir = os.environ['STUDIES']+study_id
    if not os.path.exists(study_dir):
        print 'Create directories for study '+study_id
        os.mkdir(study_dir)
        os.mkdir(study_dir+'/GeometricModel'+study_id)
        os.mkdir(study_dir+'/LVMechanics'+study_id)
        os.mkdir(study_dir+'/LVMechanics'+study_id+'/PassiveMechanics')
        os.mkdir(study_dir+'/LVMechanics'+study_id+'/ActiveMechanics')
    print 'LOG: Set up directories for '+study_id
    results_dir = os.environ['RESULTS']+study_id
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
#
#=======================================================================================================================
#


def setup_get_frames(idx):
    # This function gets the specified study ID and important frame numbers.
    filename = os.environ['PARAM_ESTIMATION'] + '/StudyNames.txt'
    f = open(filename, 'r')

    study_ids = []
    study_frames = []
    study_info = f.readline()
    while len(study_info) != 0:  # Reaching the end of the file
        study_ids.append(study_info.split()[0])
        study_frames.append(study_info.split()[1:5])
        study_info = f.readline()
        num_studies = len(study_ids)

    print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    print '      The total number of studies is', num_studies
    print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    study_id = study_ids[idx]
    study_frame = study_frames[idx]

    return study_id, study_frame
#
#=======================================================================================================================
#





