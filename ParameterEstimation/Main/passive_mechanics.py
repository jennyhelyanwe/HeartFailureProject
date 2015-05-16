__author__ = 'zwan145'

import os
from material_setup import *
from mech_setup import *
from bc_setup import *


def passive_loop_index(study_frame):
    # This function generates the MRI frame indices for the passive phase (from DS to ED).
    ds, ed, es, tot = tuple(study_frame)

    idx = [0]*(int(tot)-int(ds)+2)
    n = 0
    for i in range(int(ds), int(tot)+1):
        idx[n] = i
        n += 1
    idx[n] = int(ed)
    return idx
#
#=======================================================================================================================
#


def passive_initial_solve(study_id, study_frame, pressure, node_idx):
    # This function implements the initial forward solve from DS to ED.
    # 1. Housekeeping, set up mechanics folders.
    # 2. Loop through frames from DS to ED and for each frame:
    #    a) Get pressure and displacement boundary conditions.
    #    b) Simulate inflation.
    #    c) Save solutions in ipinit format.
    ####################################################################################################################

    ## 1. Housekeeping, set up mechanics folders.
    # Get frame numbers
    ds, ed, es, tot = tuple(study_frame)

    # Set up mechanics output folders
    dir_work = os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/PassiveMechanics'
    mech_output_setup(dir_work, 0)

    ## 2. Loop through frames from DS to ED and for each frame:
    # Set up indices for looping
    idx = passive_loop_index(study_frame)

    # Save 0 BC in warm start folder
    os.chdir(dir_work)
    copy('LV_CubicPreEpiBase_TEMPLATE.ipinit', 'WarmStartSolution/CurrentInflated_'+str(ds)+'.ipinit')

    # Loop through each passive frame and solve inflation
    for i in range(0, len(idx)-1):
        ## a) Get pressure and displacement boundary conditions.
        # Pressure increment
        p_increm = pressure[idx[i+1]-1]-pressure[idx[i]-1]
        print pressure[idx[i+1]-1], pressure[idx[i]-1]
        print 'LOG: Frame ', idx[i],' p= ',pressure[idx[i]-1]
        print 'LOG: Frame ', idx[i+1],' p= ',pressure[idx[i+1]-1]
        print 'LOG: P increment= ', p_increm

        # Update pressure BC
        bc_pressure_set(p_increm, 'WarmStartSolution/CurrentInflated_'+str(idx[i])+'.ipinit','LV_CubicPreEpiBase.ipinit')

        # Update displacement BC
        data_cur_epi = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Epi_'+str(idx[i])+'.ipdata'
        data_next_epi = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Epi_'+str(idx[i+1])+'.ipdata'
        data_cur_endo = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Endo_'+str(idx[i])+'.ipdata'
        data_next_endo = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Endo_'+str(idx[i+1])+'.ipdata'
        bc_displacement_set(node_idx, data_cur_epi, data_next_epi, data_cur_endo, data_next_endo,'LV_CubicPreEpiBase.ipinit')

        # Get geometric data
        data_epi = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Epi_'+str(idx[i+1])+'.ipdata'
        data_endo = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Endo_'+str(idx[i+1])+'.ipdata'
        copy(data_epi, 'Surface_Points_Epi_ED.ipdata')
        copy(data_endo, 'Surface_Points_Endo_ED.ipdata')

        ## b) Simulate inflation.
        os.system('cm SolveInitialInflation.com')

        ## c) Save solutions in ipinit format.
        passive_save_solutions(str(idx[i+1]), 'ForwardSolveExfile', 'ForwardSolveError', 'ForwardSolveCavityVolume',
                               'ForwardSolveStressStrain',0, False)

    # Save forward solve material parameters
    copy('LV_CubicGuc.ipmate', 'ForwardSolve.ipmate')
#
#=======================================================================================================================
#


def passive_warm_solve(study_id, frame_num):
    # This function implements the warm solve for specified frame.
    # 1. Copy current warm start solution (ipinit file) to generic name.
    # 2. Copy current frame surface data to generic name.
    # 3. Update model prediction.
    # 4. Save updated model prediction/solutions.
    ####################################################################################################################

    ## 1. Copy current warm start solution (ipinit file) to generic name.
    os.chdir(os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/PassiveMechanics')
    # Get current warm start solution
    copy('WarmStartSolution/CurrentInflated_'+frame_num+'.ipinit', 'CurrentInflated.ipinit')
    print 'WarmStartSolution/CurrentInflated_'+frame_num+'.ipinit'

    ## 2. Copy current frame surface data to generic name.
    data_epi = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Epi_'+frame_num+'.ipdata'
    data_endo = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Endo_'+frame_num+'.ipdata'
    copy(data_epi, 'Surface_Points_Epi_ED.ipdata')
    copy(data_endo, 'Surface_Points_Endo_ED.ipdata')

    ## 3. Update model prediction.
    os.system('cm SolveWarmStartPassive.com')

    ## 4. Save updated model prediction/solutions.
    passive_save_solutions(frame_num, 'OptimisedExfile', 'OptimisedError', 'OptimisedCavityVolume',
                           'OptimisedStressStrain', 1, True)
#
#=======================================================================================================================
#


def passive_save_solutions(frame_num, ex_folder, error_folder, cavity_folder, stress_folder, toggle, reg):
    # This function handles the file copying to save the results of simulations.
    # 1. Save current solution in exnode and exelem format for visualisation in CMGUI.
    # 2. Save fitting error projections for current frame for visulisation in CMGUI.
    # 3. Save stresses and strains evaluated at gauss points for current solution for visulisation in CMGUI.
    # 4. Save pressure applied at current frame for debugging checks.
    # 5. Save warm-start model prediction/solutions.
    # 5a. Save registered surface data for visualisation in CMGUI.
    # 6. Save current cavity volumes.
    ####################################################################################################################
    ## 1. Save current solution in exnode and exelem format for visualisation in CMGUI.
    copy('output/LVInflation.exnode', ex_folder+'/LVInflation_'+frame_num+'.exnode')
    copy('output/LVInflation.exelem', ex_folder+'/LVInflation_'+frame_num+'.exelem')

    ## 2. Save fitting error projections for current frame for visulisation in CMGUI.
    copy('output_errors/EndoProjectionToED.opdata', error_folder+'/EndoError_'+frame_num+'.opdata')
    copy('output_errors/EpiProjectionToED.opdata', error_folder+'/EpiError_'+frame_num+'.opdata')
    copy('output_errors/EndoProjectionToED.exdata', error_folder+'/EndoError_'+frame_num+'.exdata')
    copy('output_errors/EpiProjectionToED.exdata', error_folder+'/EpiError_'+frame_num+'.exdata')
    copy('output_errors/ED_Epi.exdata', error_folder+'/Epi_'+frame_num+'.exdata')
    copy('output_errors/ED_Endo.exdata', error_folder+'/Endo_'+frame_num+'.exdata')

    ## 3. Save stresses and strains evaluated at gauss points for current solution for visulisation in CMGUI.
    copy('output/LVInflation_gauss_ER.exdata', stress_folder+'/ER_ '+frame_num+'.exdata')
    copy('output/LVInflation_gauss_strain.exdata', stress_folder+'/Strain_'+frame_num+'.exdata')
    copy('output/LVInflation_gauss_stress.exdata', stress_folder+'/TotalStress_'+frame_num+'.exdata')
    copy('output/LVInflation_passive_gauss_stress.exdata', stress_folder+'/PassiveStress_'+frame_num+'.exdata')
    copy('output/gauss_stress.opstre', stress_folder+'/total_stress_'+frame_num+'.opstre')
    copy('output/passive_gauss_stress.opstre', stress_folder+'/passive_stress_'+frame_num+'.opstre')
    copy('output/gauss_strain.opstra', stress_folder+'/strain_'+frame_num+'.opstra')

    ## 4. Save pressure applied at current frame for debugging checks.
    copy('pressure/pressure.opvari', 'pressure/pressure_'+frame_num+'.opvari')

    if toggle == 1:
        ## 5. Save warm-start model prediction/solutions.
        copy('CurrentInflated.ipinit', 'WarmStartSolution/CurrentInflated_'+frame_num+'.ipinit')
        if reg:
            ## 5a. Save registered surface data for visualisation in CMGUI.
            copy('output_debug/ED_Endo_reg.exdata', 'OptimisedExfile/Surface_endo_reg_'+frame_num+'.exdata')
            copy('output_debug/ED_Epi_reg.exdata', 'OptimisedExfile/Surface_epi_reg_'+frame_num+'.exdata')
            copy('output_debug/ED_Endo_nonreg.exdata', 'OptimisedExfile/Surface_endo_nonreg_'+frame_num+'.exdata')
            copy('output_debug/ED_Epi_nonreg.exdata', 'OptimisedExfile/Surface_epi_nonreg_'+frame_num+'.exdata')
        ## 6. Save current cavity volumes.
        copy('output_cavity_volume/LVCavityUpdate.opelem', cavity_folder+'/LVCavity_'+frame_num+'.opelem')
    else:
        ## 5. Save warm-start model prediction/solutions.
        copy('output/LVInflation.ipinit', 'WarmStartSolution/CurrentInflated_'+frame_num+'.ipinit')
        ## 6. Save current cavity volumes.
        copy('output_cavity_volume/LVCavityCurrent.opelem', cavity_folder+'/LVCavity_'+frame_num+'.opelem')
#
#=======================================================================================================================
#
