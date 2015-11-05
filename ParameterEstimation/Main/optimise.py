__author__ = 'zwan145'

import os, sys
from passive_mechanics import *
from scipy.optimize import leastsq, fmin, fminbound, fmin_slsqp, fmin_l_bfgs_b
from scipy import array

from active_mechanics import *
from bc_setup import *
from activation_setup import *
from results import *
from datetime import datetime
import time
import math

# This python script stores functions which estimate the passive and active parameters.


def optimise_passive_main(study_id, study_frame, pressure, node_idx, forward_solve_toggle):
    # This function estimates a bulk myocardial stiffness parameter (C1) by fitting LV model predictions to
    # subject-specific image-derived geometries.
    # 1. Initial forward solve, simulating inflation from DS to ED. The model predictions at each frame is saved in
    # ipinit format.
    # 2. Generate results text file after initial solve.
    # 3. Optimise C1 to minimise summed fitting error for all frames from DS to ED.
    # 4. Evaluate final optimised solution and fitting errors.
    ####################################################################################################################

    os.chdir(os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/PassiveMechanics')
    ## 0. Customise visualisation com files.
    passive_visualisation_customise(study_id, study_frame)

    ## 0. Calculate displacements using projection of surface data onto DS.
    passive_displacement_mse(study_id, study_frame)

    ## 1. Initial forward solve, simulating inflation from DS to ED. The model predictions at each frame is saved in
    # ipinit format.
    if forward_solve_toggle == 1:
        # Initialise passive parameter
        material_create_ipmate(4.0, 'LV_CubicGuc_TEMPLATE.ipmate', 'LV_CubicGuc.ipmate')
        # Do initial solve
        passive_initial_solve(study_id, study_frame, pressure, node_idx)
        ## 2. Generate results files containing pressure, volume, endocardial and epicardial fitting errors at each frame
        # from DS+1 to ED.
        results_passive_generate(study_id, study_frame, 'ForwardSolve', 'ForwardSolve')

    ## 3. Optimise C1 to minimise summed fitting error for all frames from DS to ED.
    # Control bounds for C1 parameter.
    l_bound = 0.6
    u_bound = 20.0
    # Set initial guess for C1.
    C1_init = [4.0]
    # Copy over forward solve solutions to current warm solve solution folder.
    os.system('cp ForwardSolveSolution/*.* WarmStartSolution/.')
    # Use built-in python optimiser to estimate C1.
    C1_opt = fmin_l_bfgs_b(optimise_passive_obj_function, C1_init, approx_grad=1, bounds=[(l_bound, u_bound)],
                           epsilon=1e-2, args=[study_id, study_frame], factr=1e12)
    print C1_opt
    # Final solve using optimised C1 value.
    C1 = float(C1_opt[0][0])
    mse = optimise_passive_obj_function(C1, study_id, study_frame)

    results_passive_generate(study_id, study_frame, 'Optimised', 'Optimised')

    # Output to log file.
    print 'LOG: Final Optimised Passive Parameter: ' + str(C1)
    print 'LOG: With total MSE of '+str(mse)

    print 'LOG: Finished passive analysis. '

#
#=======================================================================================================================
#

def optimise_active_main(study_id, study_frame, pressure, node_idx):
    # This function estimates the transient of the bulk myocardial activation parameter TCa by estimating one TCa
    # parameter per MRI frame.
    # 1. Housekeeping, set up mechanics folders.
    # 2. Get ED solution from passive parameter estimation as starting point.
    # 3. Loop through all frames from ED to DS (through systole and rapid filling diastole). For each frame:
    #    a) Initial simulation using sensible guess for TCa parameter, saved model predictions in ipinit format.
    #       Convergence checking can be toggled.
    #    b) Set sensible upper and lower bounds for TCa parameter according to pressure change.
    #    c) Optimise TCa parameter for current frame to minimise fitting error at current frame.
    #    d) Evaluate final optimised solution and fitting errors.
    # 4. Generate results text file.
    ####################################################################################################################

    ## 0. Calculate displacements using projection of surface data onto DS model.
    active_displacement_mse(study_id, study_frame)

    ## 1. Housekeeping, set up mechanics folders.
    # Get important frame numbers
    ds, ed, es, tot = tuple(study_frame)

    # Set up folders
    dr = os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/ActiveMechanics'
    mech_output_setup(dr, 1) # Option 1 indicates active mechanics.

    ## 2. Get ED solution from passive parameter estimation as starting point.
    active_get_ED(study_id)

    ## 3. Loop through all frames from ED to DS (through systole and rapid filling diastole).
    # Get indices based on frame numbers
    idx = active_loop_index(study_frame)
    print idx

    # Initialise TCa
    os.chdir(dr)

    copy('TCa_ZERO_exclude_apex.ipacti', 'ForwardSolveActivation/TCa_'+str(ed)+'.ipacti')  # Zero activation at ED.
    copy('TCa_ZERO_exclude_apex.ipacti', 'OptimisedActivation/TCa_'+str(ed)+'.ipacti')  # Zero activation at ED.

    for i in range(0, len(idx)-1):
        ## a) Initial simulation using sensible guess for TCa parameter, saved model predictions in ipinit format.
        copy('WarmStartSolution/CurrentContracted_'+str(idx[i])+'.ipinit', 'CurrentContracted.ipinit')

        # Pressure increment
        p_increm = pressure[idx[i+1]-1]-pressure[idx[i]-1]
        print pressure[idx[i+1]-1], pressure[idx[i]-1]
        print 'LOG: Frame ', idx[i], ' p= ', pressure[idx[i]-1]
        print 'LOG: Frame ', idx[i+1], ' p= ', pressure[idx[i+1]-1]
        print 'LOG: P increment= ', p_increm

        # Update pressure BC
        bc_pressure_set(p_increm, 'WarmStartSolution/CurrentContracted_'+str(idx[i])+'.ipinit', 'LV_CubicPreEpiBase.ipinit')

        # Update displacement BC
        data_cur_epi = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Epi_'+str(idx[i])+'.ipdata'
        data_next_epi = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Epi_'+str(idx[i+1])+'.ipdata'
        data_cur_endo = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Endo_'+str(idx[i])+'.ipdata'
        data_next_endo = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Endo_'+str(idx[i+1])+'.ipdata'
        filename = 'OptimisedExfile/LVContraction_'+str(idx[i])+'.exnode'
        bc_displacement_set(node_idx, data_cur_epi, data_next_epi, data_cur_endo, data_next_endo,
                            'LV_CubicPreEpiBase.ipinit', filename)

        # Get geometric data
        data_epi = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Epi_'+str(idx[i+1])+'.ipdata'
        data_endo = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Endo_'+str(idx[i+1])+'.ipdata'
        copy(data_epi, 'Surface_Points_Epi_ES.ipdata')
        copy(data_endo, 'Surface_Points_Endo_ES.ipdata')

        # Get initial TCa guess
        TCa = active_init_TCa(pressure, idx[i], idx[i+1])

        # Create new ipacti file
        activation_create_ipacti_exclude_apex(TCa, 'TCa_ZERO_exclude_apex.ipacti', 'ForwardSolveActivation/TCa_'+str(idx[i+1])+'.ipacti')

        copy('ForwardSolveActivation/TCa_'+str(idx[i+1])+'.ipacti', 'TCa_current.ipacti')
        copy('OptimisedActivation/TCa_'+str(idx[i])+'.ipacti', 'TCa_previous.ipacti')

        # Monitor CMISS solution.
        CMISS_monitor = 1  # Turn CMISS process monitoring on(1)/off(0).
        if CMISS_monitor == 1:
            # Monitor the progress of the forward solve in CMISS. If program is running for too long, kill program
            # and refine load step. If program finishes, check that convergence has been reached, if not, refine load
            # step.
            converged = 0  # Default is unconverged.
            finished = 0  # Default is unfinished.
            TCa_step = 4.0  # Default TCa step is 4.0 kPa.

            # Get TCa to start from.
            with open('TCa_previous.ipacti', 'r') as file:
                data = file.readlines()
            if data[9].split()[1] == 'non-dimensional':
                TCa_prev = float(data[17].split()[-1])
            else:
                TCa_prev = float(data[21].split()[4])

            TCa_increm = TCa - TCa_prev  # Increment to be applied for current simulation.
            # Calculate number of TCa load steps.
            num_steps = abs(math.ceil(TCa_increm/TCa_step))
            if num_steps == 0:
                num_steps = 1
            active_set_load_step(TCa_step)  # Set load step by re-writing CMISS command file.

            # Set time limit for each load step increment to finish solving.
            time_limit = 160.0  # seconds.

            while converged == 0 and finished == 0:
                # Run forward solve command file in forked process.
                os.system('cm SolveInitialSystoleMF_exclude_apex.com ->& ForwardSolve.log &')
                print 'LOG: Solve to frame %d, waiting for forward solve to finish...\n' % (idx[i+1])

                # Initialise current iteration as zero.
                cur_it = 0
                with open('loaditeration.txt', 'w') as f_w:
                    f_w.write('0\n')

                # While simulation hasn't finished all load step...
                while cur_it < num_steps:
                    # Count down the time taken for one load step to finish solving. .
                    prev_it = cur_it  # Initialise previous time.
                    t_start = time.time()  # Get start time.
                    t_elapsed = time.time() - t_start  # Initialise time elapsed.
                    while t_elapsed < time_limit:
                        t_elapsed = time.time() - t_start  # Get current time elapsed.
                        count_down = time_limit - t_elapsed  # Get time remaining to finish this load step.
                        # Set up commandline timer.
                        print 'LOG: Time left: %d s\r' % count_down,
                        sys.stdout.flush()

                        # Read current load iteration number.
                        with open('loaditeration.txt', 'r') as f_i:
                            temp = f_i.readlines()
                        if temp != ['\n']:
                            if temp != []:
                                cur_it = int(temp[0])
                                # Check if CMISS has proceeded onto next load step.
                                if cur_it - prev_it > 0:
                                    print 'LOG: Iteration '+str(cur_it)+'/'+str(num_steps)+' is complete'
                                    break
                                elif cur_it == num_steps:
                                    break  # Break from while loop once all load step are complete.

                    # Once timer has run out, get current completed load step number.
                    with open('loaditeration.txt', 'r') as f_i:
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
                with open('convergence.txt', 'r') as f_c:
                    temp_c = f_c.readlines()
                    print temp_c
                if temp_c != []:
                    converged = int(temp_c[0])
                # Refine load step if did not converge.
                if converged == 0:
                    TCa_step = TCa_step * 0.5
                    num_steps = abs(math.ceil(TCa_increm/TCa_step))
                    print 'LOG: Attempt to solve contraction did not converge, refine TCa step size to '+str(TCa_step)+'\n'
                    active_set_load_step(TCa_step)
        elif CMISS_monitor == 0:
            # Run forward solve with no convergence monitoring.
            os.system('cm SolveInitialSystoleMF_exclude_apex.com')

        print 'LOG: Forward solution has finished and is converged. \n'
        time.sleep(5)
        # Save files
        active_save_solutions(str(idx[i+1]), 'ForwardSolveExfile', 'ForwardSolveError', 'ForwardSolveCavityVolume',
                              'ForwardSolveStressStrain', 0)

        ## b) Set sensible upper and lower bounds for TCa parameter according to pressure change.
        # Get bounds for TCa
        [lb, ub] = active_bounds_get(pressure, idx[i], idx[i+1])

        ## c) Optimise TCa parameter for current frame to minimise fitting error at current frame.
        TCa_Opt = fmin_l_bfgs_b(optimise_active_obj_function, [TCa], args=[study_id, str(idx[i+1])], approx_grad=1,
                                bounds=[(lb, ub)], factr=1e12, epsilon=0.1)

        print 'TCa_Opt: ' + str(TCa_Opt)

        ## d) Evaluate final optimised solution and fitting errors.
        TCa = float(TCa_Opt[0][0])
        mse = optimise_active_obj_function(TCa, study_id, str(idx[i+1]))

        print 'LOG: Optimised Passive Parameter: ' + str(TCa) + ' for frame '+str(idx[i+1])
        print 'LOG: With total MSE of ' + str(mse)

    ## 4. Generate results text file.
    results_active_generate(study_id, study_frame, 'Optimised', 'Optimised')

    ## 5. Evaluate TCa sensitivity at each optimised frame.
    active_sensitivity_evaluate(study_id, study_frame)

    print 'LOG: Finished active analysis. '
#
#=======================================================================================================================
#

def optimise_passive_obj_function(C1, study_id, study_frame):
    # This function is the objective function which interfaces with the python optimiser for passive parameter
    # estimation. It updates the C1 parameter and evaluates the overall fitting error objective function at all frames
    # between DS to ED.
    # 1. Update current C1 estimate in ipmate file.
    # 2. Re-run simulation using current C1 estimate and get new model predictions for each frame from DS+1 to ED.
    # 3. Evaluate MSE of fitting and sum for all frames from DS+1 to ED.
    # 4. Generate results text files.
    ####################################################################################################################

    print 'LOG: Evaluating C1 = '+str(C1)+'\n'

    ## 1. Update current C1 estimate in ipmate file.
    copy('LV_CubicGuc.ipmate', 'LV_CubicGuc_previous.ipmate')
    material_create_ipmate(C1, 'LV_CubicGuc_TEMPLATE.ipmate', 'LV_CubicGuc.ipmate')

    ## 2. Re-run simulation using current C1 estimate and get new model predictions for each frame from DS+1 to ED.
    idx = passive_loop_index(study_frame)
    for i in range(1, len(idx)):
        # Solve warm start using new C1 guess - using non-registration of surface data.
        passive_warm_solve(study_id, str(idx[i]))

    ## 3. Evaluate MSE of fitting and sum for all frames from DS+1 to ED.
    toggle = 1  # Scalar objective value.
    mse = optimise_passive_obj_evaluate(idx, toggle)

    print 'LOG: Finished updating solutions with current C1 = '+str(C1)
    return mse

#
#=======================================================================================================================
#

def optimise_active_obj_function(TCa, study_id, frame_num):
    # This function is the objective function which interfaces with the python optimiser for active parameter
    # estimation. It updates the TCa parameter and evaluates the fitting error objective function at the specified frame
    # 1. Update current TCa estimate.
    # 2. Re-run simulation using current TCa estimate to get new model predictions for specified frame.
    # 3. Evaluate MSE of fitting at specified frame.
    ####################################################################################################################

    ## 1. Update current TCa estimate.
    print '\033[0;30;43m LOG: Optimising... Evaluate TCa = '+str(TCa)+'\033[0m\n'
    # Update TCa parameter
    copy('OptimisedActivation/TCa_'+str(frame_num)+'.ipacti', 'TCa_previous.ipacti')
    activation_create_ipacti_exclude_apex(TCa, 'TCa_ZERO_exclude_apex.ipacti', 'TCa_current.ipacti')

    ## 2. Re-run simulation using current TCa estimate to get new model predictions for specified frame.
    active_warm_solve(study_id, frame_num)
    
    ## 3. Evaluate MSE of fitting at specified frame.
    mse = optimise_active_obj_evaluate(frame_num)

    print 'LOG: Finished updating solutions with current TCa = ' + str(TCa)
    return mse

#
#=======================================================================================================================
#

def optimise_passive_obj_evaluate(idx, toggle):
    # This function extracts the MSE of fitting for both endocardial and epicardial surfaces for all frames from DS+1 to
    # ED and sums them up.
    mse = []
    for i in range(1, len(idx)):
        # Vector error value
        info = open('OptimisedError/EpiError_'+str(idx[i])+'.opdata').read()
        array = re.split(r'[\n\\]', info)
        temp = array[3].split()
        num_epi = float(temp[len(temp)-1])
        temp = array[7].split()
        epi_rmse = float(temp[len(temp)-1])

        info = open('OptimisedError/EndoError_'+str(idx[i])+'.opdata').read()
        array = re.split(r'[\n\\]', info)
        temp = array[3].split()
        num_endo = float(temp[len(temp)-1])
        temp = array[7].split()
        endo_rmse = float(temp[len(temp)-1])

        mse.append((epi_rmse**2*num_epi + endo_rmse**2*num_endo)/(num_endo+num_epi))
        print '\033[0;30;45m LOG: Current MSE for frame '+str(idx[i])+' = '+str(mse[i-1])+'\033[0m\n'

    if toggle == 1:
        mse_tot = 0
        for i in range(1, len(idx)):
            mse_tot = mse_tot + mse[i-1]
        print 'LOG: Current total MSE for diastole = '+str(mse_tot)+'\n'
        return mse_tot
    elif toggle == 2:
        mse = array(mse)
        return mse

#
#=======================================================================================================================
#

def optimise_active_obj_evaluate(frame_num):
    # This function extracts MSE of fitting for both endocardial and epicardial surfaces and returns MSE for specified
    # frame.

    info = open('OptimisedError/EpiError_'+frame_num+'.opdata').read()
    array = re.split(r'[\n\\]', info)
    temp = array[3].split()
    num_epi = float(temp[len(temp)-1])
    temp = array[7].split()
    epi_rmse = float(temp[len(temp)-1])

    info = open('OptimisedError/EndoError_'+frame_num+'.opdata').read()
    array = re.split(r'[\n\\]', info)
    temp = array[3].split()
    num_endo = float(temp[len(temp)-1])
    temp = array[7].split()
    endo_rmse = float(temp[len(temp)-1])

    mse= (epi_rmse**2*num_epi + endo_rmse**2*num_endo)/(num_endo+num_epi)
    print '\033[0;30;45m LOG: Current MSE for frame '+frame_num+' = '+str(mse)+'\033[0m\n'

    return mse
#
#=======================================================================================================================
#

