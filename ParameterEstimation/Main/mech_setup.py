__author__ = 'zwan145'

import os
from shutil import copy

# This python script contains functions which set up the mechanical simulation.


def mech_template_setup(study_id):
    # This function copies the required mechanical simulation template CMISS files into the mechanical simulation
    # folders.
    all_f = os.listdir(os.environ['MECH_TEMPLATE']+'/PassiveMechanics')
    os.chdir(os.environ['MECH_TEMPLATE']+'/PassiveMechanics')
    dir_work = os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/PassiveMechanics/'
    for f in all_f:
        copy(f, dir_work)

    all_f = os.listdir(os.environ['MECH_TEMPLATE']+'/ActiveMechanics')
    os.chdir(os.environ['MECH_TEMPLATE']+'/ActiveMechanics')
    dir_work = os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/ActiveMechanics/'
    for f in all_f:
        copy(f, dir_work)


def mech_output_setup(dr, active_toggle):
    # This function sets up the directories within the mechanics stimulation folder which stores various outputs of
    # simulation.
    os.chdir(dr)

    if not os.path.exists('WarmStartSolution'):
        os.mkdir('WarmStartSolution')  # Store ipinit solution files.
    if not os.path.exists('ForwardSolveSolution'):
        os.mkdir('ForwardSolveSolution')  # Store forward solve ipinit solution files.

    # Folders to store outputs after initial forward solve.
    if not os.path.exists('ForwardSolveExfile'):
        os.mkdir('ForwardSolveExfile')  # Store output from inflation for visualisation.
    if not os.path.exists('ForwardSolveCavityVolume'):
        os.mkdir('ForwardSolveCavityVolume')  # Store cavity volumes after simulation.
    if not os.path.exists('ForwardSolveError'):
        os.mkdir('ForwardSolveError')  # Store error projection and RMSE.
    if not os.path.exists('ForwardSolveStressStrain'):
        os.mkdir('ForwardSolveStressStrain')  # Store stress and strain evaluations at gauss points.

    # Folders to store outputs during optimisation. Final optimised files will be here after the estimation is complete.
    if not os.path.exists('OptimisedExfile'):
        os.mkdir('OptimisedExfile')
    if not os.path.exists('OptimisedCavityVolume'):
        os.mkdir('OptimisedCavityVolume')
    if not os.path.exists('OptimisedError'):
        os.mkdir('OptimisedError')
    if not os.path.exists('OptimisedStressStrain'):
        os.mkdir('OptimisedStressStrain')

    # Active mechanical simulation requires two additional folders.
    if active_toggle == 1:
        if not os.path.exists('ForwardSolveActivation'):
            os.mkdir('ForwardSolveActivation')  # Stores ipacti file for initial forward solve.
        if not os.path.exists('OptimisedActivation'):
            os.mkdir('OptimisedActivation')  # Stores ipacti file during optimisation. Final optimised file is here.
