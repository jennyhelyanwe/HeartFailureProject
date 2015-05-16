__author__ = 'zwan145'

import os
from shutil import copy
import math
import time
import sys


def active_loop_index(study_frame):

    ds, ed, es, tot = tuple(study_frame)

    idx = [0]*int(ds)
    n = 0
    for i in range(0, int(ds)):
        idx[n] = i+1
        n += 1
    return idx
#
#=======================================================================================================================
#


def active_get_ED(study_id):
    # This function generates the MRI frame indices for the active phase (from ED to DS).
    # Get optimised files from passive mechanics folder.
    os.chdir(os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/PassiveMechanics')
    copy('LV_CubicGuc.ipmate', '../ActiveMechanics/LV_CubicGuc.ipmate')
    # Copy over ED ipinit solution.
    copy('WarmStartSolution/CurrentInflated_1.ipinit', '../ActiveMechanics/WarmStartSolution/CurrentContracted_1.ipinit')
#
#=======================================================================================================================
#


def active_init_TCa(pressure, prev_idx, cur_idx):
    # This function makes an intelligent guess about the initial TCa at current frame by considering how the pressure is
    # changing.

    # The scale by which pressure is changing is calculated.
    scale = pressure[cur_idx-1]/pressure[prev_idx-1]
    if scale == 0:
        scale = 0.5

    # Get previous TCa value.
    f = open('OptimisedActivation/TCa_'+str(prev_idx)+'.ipacti', 'r')
    data = f.readlines()
    if data[9].split()[1] == 'non-dimensional':
        TCa_prev = data[15].split()[6]
    else:
        TCa_prev = data[21].split()[4]
    f.close()

    # The initial TCa to solve for is calculated using the pressure scaling.
    if float(TCa_prev) == 0.0:
        TCa = (pressure[cur_idx-1]-pressure[prev_idx-1])*3.0
    elif (cur_idx == 3)|(cur_idx == 4):
        TCa = float(TCa_prev)*scale
    else:
        TCa = float(TCa_prev)*scale

    print 'LOG: Initial guess for TCa @ frame '+str(cur_idx)+' is '+str(TCa)+'\n'

    return TCa
#
#=======================================================================================================================
#


def active_warm_solve(study_id, frame_num):
    # This function implements the active warm solve for specified frame.
    # 1. Copy current warm start solution (ipinit file) to generic name.
    # 2. Copy current frame surface data to generic name.
    # 3. Update model prediction. Convergence can be turned on or off.
    # 4. Save updated model prediction/solutions.
    ####################################################################################################################

    ## 1. Copy current warm start solution (ipinit file) to generic name.
    os.chdir(os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/ActiveMechanics')
    copy('WarmStartSolution/CurrentContracted_'+frame_num+'.ipinit', 'CurrentContracted.ipinit')

    ## 2. Copy current frame surface data to generic name.
    data_epi = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Epi_'+frame_num+'.ipdata'
    data_endo = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Endo_'+frame_num+'.ipdata'
    copy(data_epi, 'Surface_Points_Epi_ES.ipdata')
    copy(data_endo, 'Surface_Points_Endo_ES.ipdata')

    ## 3. Update model prediction. Convergence can be turned on or off.
    # Monitor CMISS solution.
    CMISS_monitor = 1  # Turn CMISS process monitoring on(1)/off(0).
    if CMISS_monitor==1:
        # Monitor the progress of the update solve in CMISS. If program is running for too long, kill program
        # and refine load step. If program finishes, check that convergence has been reached, if not, refine load
        # step.
        converged =0  #Default is unconverged.
        TCa_step = 4.0  # Default TCa step is 4.0 kPa.


        # Get TCa to start from.
        with open('TCa_previous.ipacti', 'r') as file:
            data = file.readlines()
        if data[9].split()[1] == 'non-dimensional':
            TCa_prev = float(data[15].split()[6])
        else:
            TCa_prev = float(data[21].split()[4])

        # Get TCa to update to.
        with open('TCa_current.ipacti', 'r') as file:
            data = file.readlines()
        if data[9].split()[1] == 'non-dimensional':
            TCa = float(data[15].split()[6])
        else:
            TCa = float(data[21].split()[4])

        TCa_increm = TCa - TCa_prev  # Increment to apply for current update.
        # Calculate number of TCa load steps.
        num_steps = abs(math.ceil(TCa_increm/TCa_step))
        if num_steps == 0:
            num_steps = 1

        active_set_warm_load_step(TCa_step)  # Set load step size by re-writing CMISS command file.
        print 'LOG: Number of steps for update ' + str(num_steps) + '\n'
        # Set time limit for each load step to finish solving.
        time_limit = 60  # seconds.

        while converged == 0:
            # Run update solve command file in forked process.
            os.system('cm SolveWarmStartActive.com ->& WarmSolve.log &')
            print 'LOG: Current slide %s, waiting for TCa update to finish...\n' % frame_num

            # Initialise current iteration as zero.
            cur_it = 0
            num_steps = abs(math.ceil(TCa_increm/TCa_step))
            if num_steps == 0.0:
                num_steps = 1
            with open('warmloadstep.txt', 'w') as f_w:
                f_w.write('0\n')

            while cur_it < num_steps:
                # Count down the time taken for one load step.
                prev_it = cur_it  # Initialise previous time.
                t_start = time.time()  # Get start time.
                t_elapsed = time.time() - t_start  # Initialise time elapsed.
                while t_elapsed < time_limit:
                    t_elapsed = time.time() - t_start  # Get current time elapsed.
                    count_down = time_limit - t_elapsed  # Get time remaining to finish this step.

                    # Set up command line timer.
                    print 'LOG: Time left: %d s\r' % count_down,
                    sys.stdout.flush()

                    # Read current load iteration number.
                    with open('warmloadstep.txt', 'r') as f_i:
                        temp = f_i.readlines()
                    if temp != ['\n']:
                        if temp != []:
                            cur_it= int(temp[0])
                            # Check if CMISS has proceeded onto next load step.
                            if cur_it - prev_it > 0:
                                print 'LOG: Iteration '+str(cur_it)+'/'+str(num_steps)+' is complete'
                                break
                            elif cur_it == num_steps:
                                break  # Break from while loop once all load steps are complete.
                # Once timer has run out, get current completed load step number.
                with open('warmloadstep.txt', 'r') as f_i:
                    temp = f_i.readlines()
                if temp != []:
                    cur_it= int(temp[0])
                    # Check if current completed load step number is the final load step.
                    if (cur_it - prev_it == 0)&(cur_it!= num_steps):
                        # If current completed load step is not the final step, assume inconvergence/crash.
                        # Kill CMISS process.
                        cmd = 'ps | grep cm | grep -v grep | awk \'{print "kill -9 " $1}\' | sh'
                        os.system(cmd)
                        print 'LOG: Time ran out at iteration %d. Killed cm process.\n' % cur_it
                        break  # Break from timer while loop.
            # Check convergence.
            print 'LOG: Checking convergence! '
            time.sleep(10)  # Pause to allow text file writing to complete.
            with open('warmconvergence.txt', 'r') as f_c:
                temp_c = f_c.readlines()
            if temp_c != []:
                if temp_c != ['\n']:
                    converged = int(temp_c[0])
            # Refine load step size if it did not converge.
            if converged == 0:
                TCa_step = TCa_step * 0.5
                if TCa_step < 1e-2:
                    print 'LOG: Something has gone wrong... TCa_step = ' +str(TCa_step) + '\n'
                    quit()
                num_steps = abs(math.ceil(TCa_increm/TCa_step))
                print 'LOG: Update TCa did not converge, refine TCa step size to '+str(TCa_step)+'\n'
                active_set_warm_load_step(TCa_step)
        print 'LOG: Update TCa has converged. \n'
    elif CMISS_monitor == 0:
        # Run update solution with no convergence monitoring.
        os.system('cm SolveWarmStartActive.com')

    ## 4. Save updated model prediction/solutions.
    active_save_solutions(frame_num, 'OptimisedExfile', 'OptimisedError', 'OptimisedCavityVolume', 'OptimisedStressStrain', 1)
#
#=======================================================================================================================
#


def active_save_solutions(frame_num, ex_folder, error_folder, cavity_folder, stress_folder, toggle):
    # This function handles the file copying to save the results of simulations.
    # 1. Save current solution in exnode and exelem format for visualisation in CMGUI.
    # 2. Save fitting error projections for current frame for visulisation in CMGUI.
    # 3. Save stresses and strains evaluated at gauss points for current solution for visulisation in CMGUI.
    # 4. Save pressure applied at current frame for debugging checks.
    # 5. Save warm-start model prediction/solutions.
    # 5a. Save registered surface data for visualisation in CMGUI.
    # 6. Save current cavity volumes.
    # 7. Save activation file.
    ####################################################################################################################
    ## 1. Save current solution in exnode and exelem format for visualisation in CMGUI.
    copy('output/LVContraction.exnode', ex_folder+'/LVContraction_'+frame_num+'.exnode')
    copy('output/LVContraction.exelem', ex_folder+'/LVContraction_'+frame_num+'.exelem')
    ## 2. Save fitting error projections for current frame for visulisation in CMGUI.
    copy('output_errors/EndoProjectionToES.opdata', error_folder+'/EndoError_'+frame_num+'.opdata')
    copy('output_errors/EpiProjectionToES.opdata', error_folder+'/EpiError_'+frame_num+'.opdata')
    copy('output_errors/EndoProjectionToES.exdata', error_folder+'/EndoError_'+frame_num+'.exdata')
    copy('output_errors/EpiProjectionToES.exdata', error_folder+'/EpiError_'+frame_num+'.exdata')
    copy('output_errors/ES_Epi.exdata', error_folder+'/Epi_'+frame_num+'.exdata')
    copy('output_errors/ES_Endo.exdata', error_folder+'/Endo_'+frame_num+'.exdata')

    ## 3. Save stresses and strains evaluated at gauss points for current solution for visulisation in CMGUI.
    copy('output/LVContraction_gauss_ER.exdata', stress_folder+'/ER_'+frame_num+'.exdata')
    copy('output/LVContraction_gauss_strain.exdata', stress_folder+'/Strain_'+frame_num+'.exdata')
    copy('output/LVContraction_gauss_stress.exdata', stress_folder+'/TotalStress_'+frame_num+'.exdata')
    copy('output/LVContraction_passive_gauss_stress.exdata', stress_folder+'/PassiveStress_'+frame_num+'.exdata')
    copy('output/gauss_stress.opstre', stress_folder+'/total_stress_'+frame_num+'.opstre')
    copy('output/passive_gauss_stress.opstre', stress_folder+'/passive_stress_'+frame_num+'.opstre')
    copy('output/active_gauss_stress.opstre', stress_folder+'/active_stress_'+frame_num+'.opstre')
    copy('output/gauss_strain.opstra', stress_folder+'/strain_'+frame_num+'.opstra')
    ## 4. Save pressure applied at current frame for debugging checks.
    copy('pressure/pressure.opvari', 'pressure/pressure_'+frame_num+'.opvari')
    ## 5. Save warm-start model prediction/solutions.
    copy('output/LVContraction.ipinit', 'WarmStartSolution/CurrentContracted_'+frame_num+'.ipinit')
    if toggle == 1:
        ## 5a. Save registered surface data for visualisation in CMGUI.
        copy('output_debug/ES_Endo_reg.exdata', 'OptimisedExfile/Surface_endo_reg_'+frame_num+'.exdata')
        copy('output_debug/ES_Epi_reg.exdata', 'OptimisedExfile/Surface_epi_reg_'+frame_num+'.exdata')
        copy('output_debug/ES_Endo_nonreg.exdata', 'OptimisedExfile/Surface_endo_nonreg_'+frame_num+'.exdata')
        copy('output_debug/ES_Epi_nonreg.exdata', 'OptimisedExfile/Surface_epi_nonreg_'+frame_num+'.exdata')
        ## 6. Save current cavity volumes.
        copy('output_cavity_volume/LVCavityUpdate.opelem', cavity_folder+'/LVCavity_'+frame_num+'.opelem')
    else:
        ## 5. Save warm-start model prediction/solutions.
        copy('output/LVContraction.ipinit', 'CurrentContracted.ipinit')
        ## 6. Save current cavity volumes.
        copy('output_cavity_volume/LVCavityCurrent.opelem', cavity_folder+'/LVCavity_'+frame_num+'.opelem')

    ## 7. Save activation file.
    copy('output_debug/TCa_current.ipacti', 'OptimisedActivation/TCa_'+frame_num+'.ipacti')
#
#=======================================================================================================================
#


def active_bounds_get(pressure,  prev_idx, cur_idx):
    # This function sets an intelligent upper and lower bound for the estimation of the TCa parameter based on what the
    # LV pressure is doing.
    p = pressure[cur_idx-1] - pressure[prev_idx-1]

    f = open('OptimisedActivation/TCa_'+str(cur_idx)+'.ipacti', 'r')

    data = f.readlines()

    if data[9].split()[1] == 'non-dimensional':
        TCa_cur = data[15].split()[6]
    else:
        TCa_cur = data[21].split()[4]
    f.close()

    f = open('OptimisedActivation/TCa_'+str(prev_idx)+'.ipacti', 'r')
    data = f.readlines()
    if data[9].split()[1] == 'non-dimensional':
        TCa_prev = data[15].split()[6]
    else:
        TCa_prev = data[21].split()[4]
    f.close()

    if p > 0:
        lb = float(TCa_cur)*0.5
    else:
        lb = float(TCa_prev)*0.2

    if cur_idx == 5:
        lb = float(TCa_prev)*0.8

    ub = float(TCa_cur)*2
    if ub < 20:
        ub = 20

    print 'LOG: The bounds for frame '+str(cur_idx)+' are '+str(lb)+' and '+str(ub)+'\n'

    return lb, ub
#
#=======================================================================================================================
#


def active_set_load_step(threshold):
    # This function sets the current load step size by re-writing the CMISS command file for initial solve TCa.
    with open('SolveIVCEjectionMF.com', 'r') as file:
        data = file.readlines()

    data[38] = '$MAXIMUM_INCREM = ceil(abs($TCa_step)/'+str(threshold)+');\n'

    with open('SolveIVCEjectionMF.com', 'w') as file:
        file.writelines(data)
#
#=======================================================================================================================
#


def active_set_warm_load_step(threshold):
    # This function set the current load step size by re-writing the CMISS command file for updating TCa.
    with open('SolveWarmStartActive.com', 'r') as file:
        data = file.readlines()

    data[42] = '$MAXIMUM_INCREM = ceil(abs($TCa_increm)/'+str(threshold)+');\n'

    with open('SolveWarmStartActive.com', 'w') as file:
        file.writelines(data)
#
#=======================================================================================================================
#
