__author__ = 'zwan145'

import os, sys
import scipy
from scipy.optimize import leastsq, fmin, fminbound, fmin_slsqp, fmin_l_bfgs_b, fmin_cobyla
from numpy import array
from results import *
import time
import math
from passive_mechanics import *
import nlopt

# Global parameters to store trajectory of optimisation.
parameter_trajectory = []
objective_function_trajectory = []

class bcolors:
    HEADER = '\033[0;30;45m'
    OKBLUE = '\033[0;30;44m'
    OKGREEN = '\033[0;30;42m'
    WARNING = '\033[0;30;43m'
    FAIL = '\033[0;30;41m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# This python script stores functions which estimate the passive and active parameters.
def optimise_passive_main(study_id, study_frame, pressure, node_idx, C1_init, forward_solve_toggle):
    # This function estimates a bulk myocardial stiffness parameter (C1) by fitting LV model predictions to
    # subject-specific image-derived geometries.
    # 1. Initial forward solve, simulating inflation from DS to ED. The model predictions at each frame is saved in
    # ipinit format.
    # 2. Generate results text file after initial solve.
    # 3. Optimise C1 to minimise summed fitting error for all frames from DS to ED.
    # 4. Evaluate final optimised solution and fitting errors.
    ####################################################################################################################
    print 'LOG_STATUS: Begin optimised_passive_main.\n'
    os.chdir(os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/PassiveMechanics')
    # 0. Customise visualisation com files.
    passive_visualisation_customise(study_id, study_frame)

    # 0. Calculate displacements using projection of surface data onto DS.
    #passive_displacement_mse(study_id, study_frame)
    # 1. Initial forward solve, simulating inflation from DS to ED. The model predictions at each frame is saved in
    # ipinit format.
    if forward_solve_toggle == 1:
        # Initialise passive parameter
        material_create_ipmate(C1_init, 'LV_CubicGuc_TEMPLATE.ipmate', 'LV_CubicGuc.ipmate')
        # Do initial solve
        passive_initial_solve(study_id, study_frame, pressure, node_idx)
        # 2. Generate results files containing pressure, volume, endocardial and epicardial fitting errors at each frame
        # from DS+1 to ED.
        results_passive_generate(study_id, study_frame, 'ForwardSolve', 'ForwardSolve')
        optimise_passive_exclude_data(study_frame)

    # 3. Optimise C1 to minimise summed fitting error for all frames from DS to ED.
    # Control bounds for C1 parameter.
    #l_bound = 0.3
    #u_bound = 30.0
    # Copy over forward solve solutions to current warm solve solution folder.
    os.system('cp ForwardSolveSolution/*.* WarmStartSolution/.')
    copy('LV_CubicGuc_ForwardSolve.ipmate', 'LV_CubicGuc.ipmate')

    # Use NLOpt COBYLA optimiser to evaluate C1.
    idx = passive_loop_index(study_frame)
    toggle = 1  # For scalar objective function.
    opt = nlopt.opt(nlopt.LN_COBYLA, 1)
    opt.set_min_objective(lambda x, grad: optimise_passive_obj_function(x, study_id, study_frame, toggle))
    opt.set_lower_bounds(0.3)
    opt.set_upper_bounds(15)
    opt.set_ftol_abs(0.01*len(idx))
    opt.set_initial_step(0.1)
    opt.set_maxeval(100)
    C1_init = [C1_init]
    c1_opt = opt.optimize(C1_init)
    opt_val = opt.last_optimum_value()
    print bcolors.OKGREEN + 'Optimiser log: The optimised C1 value is ' + str(c1_opt) + bcolors.ENDC
    print bcolors.OKGREEN + 'Optimiser log: The optimised objective function value is ' + str(opt_val) + bcolors.ENDC
    # Get optimised results

    result = opt.last_optimize_result()
    if result > 0:
        print bcolors.OKGREEN + 'Optimiser log: Successful termination' + bcolors.ENDC
        if result == 1:
            print bcolors.OKGREEN + 'Optimiser log: stopped with generic success.'+ bcolors.ENDC
        elif result == 2:
            print bcolors.OKGREEN + 'Optimiser log: stopped because stopval was reached.'+ bcolors.ENDC
        elif result == 3:
            print bcolors.OKGREEN + 'Optimiser log: stopped because ftol_rel or ftol_abs was reached.'+ bcolors.ENDC
        elif result == 4:
            print bcolors.OKGREEN + 'Optimiser log: stopped because xtol_rel or xtol_abs was reached.'+ bcolors.ENDC
        elif result == 5:
            print bcolors.WARNING + 'Optimiser log: stopped because maxeval was reached.' + bcolors.ENDC
        elif result == 6:
            print bcolors.WARNING + 'Optimiser log: stopped because maxtime was reached.' + bcolors.ENDC
    else:
        print bcolors.FAIL + 'Optimiser log: Error code termination' + bcolors.ENDC
        if result == -1:
            print bcolors.FAIL + 'Optimiser log: stopped with generic failure.' + bcolors.ENDC
        elif result == -2:
            print bcolors.FAIL + 'Optimiser log: invalid arguments - e.g. lower bounds are bigger than upper bounds, ' \
                                 'or unknown algorithm.' + bcolors.ENDC
        elif result == -3:
            print bcolors.FAIL + 'Optimiser log: ran out of memory.' + bcolors.ENDC
        elif result == -4:
            print bcolors.FAIL + 'Optimiser log: halted due to roundoff errors limiting progress. ' + bcolors.ENDC
        elif result == -5:
            print bcolors.FAIL + 'Optimiser log: halted because of forced termination - user called nlopt_force_stop(' \
                                 'opt).' + bcolors.ENDC
    mse = optimise_passive_obj_function(c1_opt, study_id, study_frame, toggle)
    if (mse - opt_val) < 1e-2:
        print 'Double checked objective function at optimised C1. '
    else:
        print bcolors.FAIL + 'ERROR: Objective function at optimal C1: '+str(mse)+' doesn''t match what optimiser ' \
                                                                                  'gives: '+str(opt_val) + bcolors.ENDC
        quit()
    results_passive_generate(study_id, study_frame, 'Optimised', 'Optimised')
    global parameter_trajectory
    global objective_function_trajectory
    with open('OptimisationTrajectory_'+str(study_id)+'.txt', 'w') as f:
        f.write('C1\tMSE\'n')
        for i in range(0, len(parameter_trajectory)):
            f.write(str(parameter_trajectory[i])+ '\t'+str(objective_function_trajectory[i])+ '\n')
    return c1_opt, opt_val


def optimise_active_main(study_id, study_frame, pressure, node_idx, apex_acti):
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
    print 'LOG_STATUS: Begin optimised_active_main.\n'
    # 0. Calculate displacements using projection of surface data onto DS model.
    active_displacement_mse(study_id, study_frame)

    # 1. Housekeeping, set up mechanics folders.
    # Get important frame numbers
    ds, ed, es, tot = tuple(study_frame)

    # Set up folders
    dr = os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/ActiveMechanics'
    mech_output_setup(dr, 1)  # Option 1 indicates active mechanics.

    # 2. Get ED solution from passive parameter estimation as starting point.
    active_get_ED(study_id)

    # 3. Loop through all frames from ED to DS (through systole and rapid filling diastole).
    # Get indices based on frame numbers
    idx = active_loop_index(study_frame)
    print idx

    # Initialise TCa
    os.chdir(dr)

    copy('TCa_ZERO_exclude_apex.ipacti', 'ForwardSolveActivation/TCa_'+str(ed)+'.ipacti')  # Zero activation at ED.
    copy('TCa_ZERO_exclude_apex.ipacti', 'OptimisedActivation/TCa_'+str(ed)+'.ipacti')  # Zero activation at ED.

    for i in range(0, len(idx)-1):
        print 'LOG_STATUS_CURRENT: Begin frame '+str(idx[i])+'\n'
        # a) Initial simulation using sensible guess for TCa parameter, saved model predictions in ipinit format.
        copy('WarmStartSolution/CurrentContracted_'+str(idx[i])+'.ipinit', 'CurrentContracted.ipinit')

        # Pressure increment
        p_increm = pressure[idx[i+1]-1]-pressure[idx[i]-1]
        print pressure[idx[i+1]-1], pressure[idx[i]-1]
        print 'LOG: Frame ', idx[i], ' p= ', pressure[idx[i]-1]
        print 'LOG: Frame ', idx[i+1], ' p= ', pressure[idx[i+1]-1]
        print 'LOG: P increment= ', p_increm

        # Update pressure BC
        bc_pressure_set(p_increm, 'WarmStartSolution/CurrentContracted_'+str(idx[i])+'.ipinit',
                        'LV_CubicPreEpiBase.ipinit')

        # Update displacement BC
        data_cur_epi = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Epi_'+str(idx[i])+'.ipdata'
        data_next_epi = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Epi_'+str(
            idx[i+1]) + '.ipdata'
        data_cur_endo = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Endo_'+str(
            idx[i]) + '.ipdata'
        data_next_endo = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Endo_'+str(
            idx[i+1]) + '.ipdata'
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
        if apex_acti == 1:
            TCa_apex = TCa
        else:
            TCa_apex = 0.0
        activation_create_ipacti(TCa, TCa_apex, 'TCa_ZERO_exclude_apex.ipacti', 'ForwardSolveActivation/TCa_'
                                 + str(idx[i+1])+'.ipacti')

        copy('ForwardSolveActivation/TCa_'+str(idx[i+1])+'.ipacti', 'TCa_current.ipacti')
        copy('OptimisedActivation/TCa_'+str(idx[i])+'.ipacti', 'TCa_previous.ipacti')

        # Monitor CMISS solution.
        CMISS_monitor = 1  # Turn CMISS process monitoring on(1)/off(0).
        if CMISS_monitor == 1:
            # Monitor the progress of the forward solve in CMISS. If program is running for too long, kill program
            # and refine load step. If program finishes, check that convergence has been reached, if not, refine load
            # step.
            converged = 0  # Default is not converged.
            finished = 0  # Default is unfinished.
            TCa_step = 4.0  # Default TCa step is 4.0 kPa.

            # Get TCa to start from.
            with open('TCa_previous.ipacti', 'r') as file_:
                data = file_.readlines()
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
                os.system('cm SolveInitialSystoleMF.com ->& ForwardSolve.log &')
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
                        if (cur_it - prev_it == 0) & (cur_it != num_steps):
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
                    TCa_step *= 0.5
                    if TCa_step < 0.1:
                        print 'LOG_STATUS_CURRENT: Something has gone wrong at (forward solve) frame '+str(idx[i])+'... TCa_step = ' \
                              + str(TCa_step) + '\n'
                        quit()
                    num_steps = abs(math.ceil(TCa_increm/TCa_step))
                    print 'LOG: Attempt to solve contraction did not converge, refine TCa step size to '+str(TCa_step)\
                          + '\n'
                    active_set_load_step(TCa_step)
        elif CMISS_monitor == 0:
            # Run forward solve with no convergence monitoring.
            os.system('cm SolveInitialSystoleMF.com')

        print 'LOG: Forward solution has finished and is converged. \n'
        time.sleep(5)
        # Save files
        active_save_solutions(str(idx[i+1]), 'ForwardSolveExfile', 'ForwardSolveError', 'ForwardSolveCavityVolume',
                              'ForwardSolveStressStrain', 0)

        # b) Set sensible upper and lower bounds for TCa parameter according to pressure change.
        # Get bounds for TCa
        [lb, ub] = active_bounds_get(pressure, idx[i], idx[i+1])

        # c) Optimise TCa parameter for current frame to minimise fitting error at current frame.
        TCa_Opt = fmin_l_bfgs_b(optimise_active_obj_function, [TCa], args=[study_id, str(idx[i+1]), apex_acti],
                                approx_grad=1, bounds=[(lb, ub)], factr=1e12, epsilon=0.1)

        print 'TCa_Opt: ' + str(TCa_Opt)

        # d) Evaluate final optimised solution and fitting errors.
        TCa = float(TCa_Opt[0][0])
        mse = optimise_active_obj_function(TCa, study_id, str(idx[i+1]), apex_acti)

        print 'LOG_STATUS: Optimised Passive Parameter: ' + str(TCa) + ' for frame '+str(idx[i+1])
        print 'LOG_STATUS: With total MSE of ' + str(mse)
        print 'LOG_STATUS_CURRENT: Completed frame '+str(idx[i])+'\n'
    # 4. Generate results text file.

    results_active_generate(study_id, study_frame, 'Optimised', 'Optimised')
    print 'LOG_STATUS: Completed optimised_active_main.\n'


def optimise_passive_obj_function(C1, study_id, study_frame, toggle):
    # This function is the objective function which interfaces with the python optimiser for passive parameter
    # estimation. It updates the C1 parameter and evaluates the overall fitting error objective function at all frames
    # between DS to ED.
    # 1. Update current C1 estimate in ipmate file.
    # 2. Re-run simulation using current C1 estimate and get new model predictions for each frame from DS+1 to ED.
    # 3. Evaluate MSE of fitting and sum for all frames from DS+1 to ED.
    # 4. Generate results text files.
    ####################################################################################################################
    os.chdir(os.environ['STUDIES'] + study_id + '/LVMechanics' + study_id + '/PassiveMechanics')
    print 'LOG_STATUS: Evaluating C1 = '+str(C1)+'\n'
    # 1. Update current C1 estimate in ipmate file.
    copy('LV_CubicGuc.ipmate', 'LV_CubicGuc_previous.ipmate')
    material_create_ipmate(C1, 'LV_CubicGuc_TEMPLATE.ipmate', 'LV_CubicGuc.ipmate')

    # 2. Re-run simulation using current C1 estimate and get new model predictions for each frame from DS+1 to ED.
    idx = passive_loop_index(study_frame)
    for i in range(1, len(idx)):
        # Solve warm start using new C1 guess - using non-registration of surface data.
        passive_warm_solve(study_id, str(idx[i]))

    # 3. Evaluate MSE of fitting and sum for all frames from DS+1 to ED.
    mse = optimise_passive_obj_evaluate_sensible_projection(idx,  toggle)
    global C1_previous
    C1_previous = C1
    global parameter_trajectory
    parameter_trajectory.append(C1)
    global objective_function_trajectory
    objective_function_trajectory.append(mse)
    print 'LOG_STATUS: Finished updating solutions with current C1 = '+str(C1)
    return mse


def optimise_passive_obj_derivative(C1, study_id, study_frame, toggle):
    # This function evaluates the gradient around the objective function using the epsilon provided and
    # using central difference approximation.
    epsilon = 1e-2
    print 'LOG_STATUS: Calculating first derivative around current C1 = '+str(C1)+' with epsilon = '+str(epsilon)+'\n'
    print 'LOG_STATUS: Evaluating C1 = '+str(C1+epsilon/2)+'\n'
    copy('LV_CubicGuc.ipmate', 'LV_CubicGuc_previous.ipmate')
    material_create_ipmate(C1+epsilon/2, 'LV_CubicGuc_TEMPLATE.ipmate', 'LV_CubicGuc.ipmate')

    # 2. Re-run simulation using current C1 estimate and get new model predictions for each frame from DS+1 to ED.
    idx = passive_loop_index(study_frame)
    for i in range(1, len(idx)):
        # Solve warm start using new C1 guess - using non-registration of surface data.
        passive_warm_solve(study_id, str(idx[i]))

    # 3. Evaluate MSE of fitting and sum for all frames from DS+1 to ED.
    mse_plus = optimise_passive_obj_evaluate(idx, toggle)

    print 'LOG_STATUS: Evaluating C1 = '+str(C1-epsilon/2)+'\n'
    copy('LV_CubicGuc.ipmate', 'LV_CubicGuc_previous.ipmate')
    material_create_ipmate(C1-epsilon/2, 'LV_CubicGuc_TEMPLATE.ipmate', 'LV_CubicGuc.ipmate')

    # 2. Re-run simulation using current C1 estimate and get new model predictions for each frame from DS+1 to ED.
    idx = passive_loop_index(study_frame)
    for i in range(1, len(idx)):
        # Solve warm start using new C1 guess - using non-registration of surface data.
        passive_warm_solve(study_id, str(idx[i]))

    # 3. Evaluate MSE of fitting and sum for all frames from DS+1 to ED.
    mse_minus = optimise_passive_obj_evaluate(idx, toggle)

    first_derivative = float((mse_plus - mse_minus)/epsilon)
    first_derivative = numpy.array([first_derivative])
    return first_derivative


def optimise_active_obj_function(TCa, study_id, frame_num, apex_acti):
    # This function is the objective function which interfaces with the python optimiser for active parameter
    # estimation. It updates the TCa parameter and evaluates the fitting error objective function at the specified frame
    # 1. Update current TCa estimate.
    # 2. Re-run simulation us ing current TCa estimate to get new model predictions for specified frame.
    # 3. Evaluate MSE of fitting at specified frame.
    ####################################################################################################################

    # 1. Update current TCa estimate.
    print '\033[0;30;43m LOG: Optimising... Evaluate TCa = '+str(TCa)+'\033[0m\n'
    # Update TCa parameter
    copy('OptimisedActivation/TCa_'+str(frame_num)+'.ipacti', 'TCa_previous.ipacti')
    if apex_acti == 1:
        TCa_apex = TCa
    else:
        TCa_apex = 0.0
    activation_create_ipacti(TCa, TCa_apex, 'TCa_ZERO_exclude_apex.ipacti', 'TCa_current.ipacti')

    # 2. Re-run simulation using current TCa estimate to get new model predictions for specified frame.
    active_warm_solve(study_id, frame_num)

    # 3. Evaluate MSE of fitting at specified frame.
    mse = optimise_active_obj_evaluate(frame_num)

    print 'LOG: Finished updating solutions with current TCa = ' + str(TCa) + ' and TCa_apex = ' + str(TCa_apex)+'\n'
    return mse


def optimise_passive_obj_evaluate(idx, toggle):
    # This function extracts the MSE of fitting for both endocardial and epicardial surfaces for all frames from DS+1 to
    # ED and sums them up.
    mse = []
    for i in range(1, len(idx)):
        # Vector error value
        info = open('OptimisedError/EpiError_'+str(idx[i])+'.opdata').read()
        array_ = re.split(r'[\n\\]', info)
        temp = array_[3].split()
        num_epi = float(temp[len(temp)-1])
        temp = array_[7].split()
        epi_rmse = float(temp[len(temp)-1])

        info = open('OptimisedError/EndoError_'+str(idx[i])+'.opdata').read()
        array_ = re.split(r'[\n\\]', info)
        temp = array_[3].split()
        num_endo = float(temp[len(temp)-1])
        temp = array_[7].split()
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


def _read_exdata_values(filename):
    # Read exdata
    data_num = []
    value = []
    coord = []
    with open(filename, 'r') as f:
        for i in range(0, 10):
            junk = f.readline()
        temp = f.readline()
        while temp != '':
            data_num.append(temp.split()[-1])
            temp = f.readline()
            coord.append([float(temp.split()[0]), float(temp.split()[1]), float(temp.split()[2])])
            temp = f.readline()
            value.append([float(temp.split()[0]), float(temp.split()[1]), float(temp.split()[2])])
            temp = f.readline()
    return coord, value


def _read_exdata_stress_values(filename):
    # Read exdata
    elem = []
    value = []
    xi = []
    with open(filename, 'r') as f:
        for i in range(0, 16):
            junk = f.readline()
        temp = f.readline()
        while temp != '':
            temp = f.readline()
            elem.append(int(temp.split()[1]))
            xi_1 = float(temp.split()[3])
            xi_2 = float(temp.split()[4])
            xi_3 = float(temp.split()[5])
            xi.append([xi_1, xi_2, xi_3])
            temp = f.readline()
            value_pt = []
            for j in range(0, 6):
                value_pt.append(float(temp.split()[0]))
                temp = f.readline()
            value.append(value_pt)
    return elem, xi, value


def _write_exdata_values(filename, coords, values, header):
    with open(filename, 'w') as f:
        f.write(' Group name: '+header+'\n')
        f.write(' #Fields=2\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 1, #Derivatives=0\n')
        f.write('   y.  Value index= 2, #Derivatives=0\n')
        f.write('   z.  Value index= 3, #Derivatives=0\n')
        f.write(' 2) error, field, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 4, #Derivatives=0\n')
        f.write('   y.  Value index= 5, #Derivatives=0\n')
        f.write('   z.  Value index= 6, #Derivatives=0\n')
        for i in range(0, len(coords)):
            f.write(' Node:\t'+str(i+1)+'\n')
            f.write(str(coords[i][0])+'\t'+str(coords[i][1])+'\t'+str(coords[i][2])+'\n')
            f.write(str(values[i][0])+'\t'+str(values[i][1])+'\t'+str(values[i][2])+'\n')


def _read_xi_values(filename):
    elem = []
    xi = []
    with open(filename, 'r') as f:
        temp = f.readline()
        while temp != '':
            elem.append(int(temp.split()[1]))
            xi.append([float(temp.split()[2]), float(temp.split()[3]), float(temp.split()[4])])
            temp = f.readline()
    return elem, xi


def optimise_passive_exclude_data(study_frame):
    idx = passive_loop_index(study_frame)
    for i in range(1, len(idx)):
        [elem_epi, xi_epi] = _read_xi_values('ForwardSolveError/ProjectionXiEpi_' + str(idx[i]) + '.ipxi')
        [elem_endo, xi_endo] = _read_xi_values('ForwardSolveError/ProjectionXiEndo_' + str(idx[i]) + '.ipxi')
        # Concatenate magnitudes of error projection epicardium, excluding those at border of basal plane.
        TOL = 1e-2
        data_points_excluded_epi = []
        for j in range(0, len(elem_epi)):
            if elem_epi[j] > 12:
                if abs(xi_epi[j][1] - 1.0) < TOL:
                    data_points_excluded_epi.append(j)
        # Concatenate magnitude of error projection for endocardium, excluding those at border of basal plane.
        data_points_excluded_endo = []
        for j in range(0, len(elem_endo)):
            if elem_endo[j] > 12:
                if abs(xi_endo[j][1] - 1.0) < TOL:
                    data_points_excluded_endo.append(j)
    numpy.savetxt('ForwardSolveError/DataPointsExcludedEpi.txt', data_points_excluded_epi)
    numpy.savetxt('ForwardSolveError/DataPointsExcludedEndo.txt', data_points_excluded_endo)


def optimise_passive_obj_evaluate_sensible_projection(idx, toggle):
    data_points_excluded_epi = numpy.loadtxt('ForwardSolveError/DataPointsExcludedEpi.txt', float)
    data_points_excluded_endo = numpy.loadtxt('ForwardSolveError/DataPointsExcludedEndo.txt', float)
    error_magnitude = []
    num_epi_excluded = []
    num_endo_excluded = []
    for i in range(1, len(idx)):
        coord_epi_updated = []
        coord_endo_updated = []
        projection_epi_updated = []
        projection_endo_updated = []
        [coord_epi, projection_epi] = _read_exdata_values('OptimisedError/EpiError_'+str(idx[i])+'.exdata')
        [coord_endo, projection_endo] = _read_exdata_values('OptimisedError/EndoError_'+str(idx[i])+'.exdata')
        # Concatenate magnitudes of error projection epicardium, excluding those at border of basal plane.
        error_magnitude_epi = []
        for j in range(0, len(coord_epi)):
            if not j in data_points_excluded_epi:
                error_magnitude_epi.append(numpy.linalg.norm(projection_epi[j][0:3])**2)
                coord_epi_updated.append(coord_epi[j][0:3])
                projection_epi_updated.append(projection_epi[j][0:3])
        # Concatenate magnitude of error proejction for endocardium, excluding those at border of basal plane.
        error_magnitude_endo = []
        for j in range(0, len(coord_endo)):
            if not j in data_points_excluded_endo:
                error_magnitude_endo.append(numpy.linalg.norm(projection_endo[j][0:3])**2)
                coord_endo_updated.append(coord_endo[j][0:3])
                projection_endo_updated.append(projection_endo[j][0:3])
        # Evaluate MSE for frame.
        num_epi_excluded.append(len(data_points_excluded_epi))
        num_endo_excluded.append(len(data_points_excluded_endo))
        error_magnitude.append((sum(error_magnitude_endo) +
                                sum(error_magnitude_epi))/(len(error_magnitude_endo)+len(error_magnitude_epi)))
        # Rewrite error projection files
        print len(coord_epi_updated)
        print len(projection_epi_updated)
        _write_exdata_values('OptimisedError/EpiError_'+str(idx[i])+'_sensible.exdata', coord_epi_updated,
                             projection_epi_updated, 'EpiProjectionToED')
        _write_exdata_values('OptimisedError/EndoError_' + str(idx[i]) + '_sensible.exdata', coord_endo_updated,
                             projection_endo_updated, 'EndoProjectionToED')
    if toggle == 1:
        total_error_magnitude = sum(error_magnitude)
        print 'LOG: Current total error of projection of all passive frames = ' + str(total_error_magnitude) + '\n'
        print 'LOG: Number of data points excluded for epicardium: '
        print num_epi_excluded
        print 'LOG: Number of data points excluded for endocardium: '
        print num_endo_excluded
        return total_error_magnitude
    elif toggle == 2:
        print 'LOG: Largest total error of projection is ' + str(max(error_magnitude)) + '\n'
        return error_magnitude


def optimise_passive_obj_evaluate_STF_10_only(idx, toggle):
    error_magnitude = []
    num_epi_excluded = []
    num_endo_excluded = []

    for i in range(1, len(idx)):
        coord_epi_updated = []
        coord_endo_updated = []
        projection_epi_updated = []
        projection_endo_updated = []
        [coord_epi, projection_epi] = _read_exdata_values('OptimisedError/EpiError_'+str(idx[i])+'.exdata')
        [coord_endo, projection_endo] = _read_exdata_values('OptimisedError/EndoError_'+str(idx[i])+'.exdata')
        [elem_epi, xi_epi] = _read_xi_values('OptimisedError/ProjectionXiEpi_'+str(idx[i])+'.ipxi')
        [elem_endo, xi_endo] = _read_xi_values('OptimisedError/ProjectionXiEndo_' + str(idx[i]) + '.ipxi')

        # Concatenate magnitudes of error projection epicardium, excluding those at border of basal plane.
        TOL = 1e-2
        error_magnitude_epi = []
        data_points_excluded_epi = []
        for j in range(0, len(elem_epi)):
            if elem_epi[j] > 12:
                if abs(xi_epi[j][1] - 1.0) > TOL:
                    # Include data point (not projecting to edge of basal element).
                    error_magnitude_epi.append(numpy.linalg.norm(projection_epi[j][0:3])**2)
                    coord_epi_updated.append(coord_epi[j][0:3])
                    projection_epi_updated.append(projection_epi[j][0:3])
                else:
                    data_points_excluded_epi.append(j)
            elif elem_epi[j] in [6, 7, 10, 11]:
                data_points_excluded_epi.append(j)
            else:
                error_magnitude_epi.append(numpy.linalg.norm(projection_epi[j][0:3])**2)
                coord_epi_updated.append(coord_epi[j][0:3])
                projection_epi_updated.append(projection_epi[j][0:3])

        # Concatenate magnitude of error proejction for endocardium, excluding those at border of basal plane.
        error_magnitude_endo = []
        data_points_excluded_endo = []
        for j in range(0, len(elem_endo)):
            if elem_endo[j] > 12:
                if abs(xi_endo[j][1] - 1.0) > TOL:
                    # Include data point (not projecting to edge of basal element).
                    error_magnitude_endo.append(numpy.linalg.norm(projection_endo[j][0:3])**2)
                    coord_endo_updated.append(coord_endo[j][0:3])
                    projection_endo_updated.append(projection_endo[j][0:3])
                else:
                    data_points_excluded_endo.append(j)
            elif elem_epi[j] in [6, 7, 10, 11]:
                data_points_excluded_endo.append(j)
            else:
                error_magnitude_endo.append(numpy.linalg.norm(projection_endo[j][0:3])**2)
                coord_endo_updated.append(coord_endo[j][0:3])
                projection_endo_updated.append(projection_endo[j][0:3])
        # Evaluate MSE for frame.
        num_epi_excluded.append(len(data_points_excluded_epi))
        num_endo_excluded.append(len(data_points_excluded_endo))
        error_magnitude.append((sum(error_magnitude_endo) +
                                sum(error_magnitude_epi))/(len(error_magnitude_endo)+len(error_magnitude_epi)))
        # Rewrite error projection files
        print len(coord_epi_updated)
        print len(projection_epi_updated)
        #print 'OptimisedError/EpiError_' + str(idx[i]) + '.exdata'
        _write_exdata_values('OptimisedError/EpiError_'+str(idx[i])+'_sensible.exdata', coord_epi_updated,
                             projection_epi_updated, 'EpiProjectionToED')
        _write_exdata_values('OptimisedError/EndoError_' + str(idx[i]) + '_sensible.exdata', coord_endo_updated,
                             projection_endo_updated, 'EndoProjectionToED')

    if toggle == 1:
        total_error_magnitude = sum(error_magnitude)
        print 'LOG: Current total error of projection of all passive frames = ' + str(total_error_magnitude) + '\n'
        print 'LOG: Number of data points excluded for epicardium: '
        print num_epi_excluded
        print 'LOG: Number of data points excluded for endocardium: '
        print num_endo_excluded
        return total_error_magnitude
    elif toggle == 2:
        print 'LOG: Largest total error of projection is ' + str(max(error_magnitude)) + '\n'
        return error_magnitude


def optimise_active_obj_evaluate(frame_num):
    # This function extracts MSE of fitting for both endocardial and epicardial surfaces and returns MSE for specified
    # frame.

    info = open('OptimisedError/EpiError_'+frame_num+'.opdata').read()
    array_ = re.split(r'[\n\\]', info)
    temp = array_[3].split()
    num_epi = float(temp[len(temp)-1])
    temp = array_[7].split()
    epi_rmse = float(temp[len(temp)-1])

    info = open('OptimisedError/EndoError_'+frame_num+'.opdata').read()
    array_ = re.split(r'[\n\\]', info)
    temp = array_[3].split()
    num_endo = float(temp[len(temp)-1])
    temp = array_[7].split()
    endo_rmse = float(temp[len(temp)-1])
    mse = (epi_rmse**2*num_epi + endo_rmse**2*num_endo)/(num_endo+num_epi)
    print '\033[0;30;45m LOG: Current MSE for frame '+frame_num+' = '+str(mse)+'\033[0m\n'
    return mse


def optimise_determine_apical_activation(study_id, study_frame, pressure, node_idx):
    print 'LOG_STATUS: Begin optimise_determine_apical_activation.\n'
    time.sleep(10)
    # This function solves the forward solution to end systole using either apical activation or no apical activation.
    # For both cases a TCa value is optimised, and the optimised MSE is compared to determine whether to use one method
    # or the other.

    ds, ed, es, tot = tuple(study_frame)
    # Take TCa initial guesses to be 80 kPa.
    active_get_ED(study_id)
    dr = os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/ActiveMechanics'
    os.chdir(dr)

    copy('TCa_ZERO_exclude_apex.ipacti', 'ForwardSolveActivation/TCa_'+str(ed)+'.ipacti')  # Zero activation at ED.
    copy('TCa_ZERO_exclude_apex.ipacti', 'OptimisedActivation/TCa_'+str(ed)+'.ipacti')  # Zero activation at ED.

    copy('WarmStartSolution/CurrentContracted_1.ipinit', 'CurrentContracted.ipinit')

    # Pressure increment
    p_increm = pressure[int(es)-1]-pressure[0]
    print pressure[int(es)-1], pressure[0]
    print 'LOG: ED Frame ', ed, ' p= ', pressure[0]
    print 'LOG: ES Frame ', es, ' p= ', pressure[int(es)-1]
    print 'LOG: P increment= ', p_increm

    # Update pressure BC
    bc_pressure_set(p_increm, 'WarmStartSolution/CurrentContracted_'+str(ed)+'.ipinit', 'LV_CubicPreEpiBase.ipinit')

    # Update displacement BC
    data_cur_epi = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Epi_'+str(ed)+'.ipdata'
    data_next_epi = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Epi_'+str(es)+'.ipdata'
    data_cur_endo = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Endo_'+str(ed)+'.ipdata'
    data_next_endo = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Endo_'+str(es)+'.ipdata'
    filename = 'OptimisedExfile/LVContraction_'+str(ed)+'.exnode'
    bc_displacement_set(node_idx, data_cur_epi, data_next_epi, data_cur_endo, data_next_endo,
                        'LV_CubicPreEpiBase.ipinit', filename)

    # Get geometric data
    data_epi = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Epi_'+str(es)+'.ipdata'
    data_endo = os.environ['GEOM_DATA']+study_id+'/Active/'+study_id+'_Surface_Points_Endo_'+str(es)+'.ipdata'
    copy(data_epi, 'Surface_Points_Epi_ES.ipdata')
    copy(data_endo, 'Surface_Points_Endo_ES.ipdata')

    # Set initial guess of TCa to be 5 times that of the es pressure.
    TCa = pressure[int(es)-1]*5.5
    TCa_apex = [TCa, 0]  # First optimisation uses apical activation, second optimisation excludes apical activation.
    mse = [0, 0]  # MSE array used to store the mse of two optimisations.
    for i in [0, 1]:
        if i == 1:
            print 'LOG_STATUS_CURRENT: Optimising with apical activation.\n'
            time.sleep(5)
        else:
            print 'LOG_STATUS_CURRENT: Optimising without apical activation. \n'
        activation_create_ipacti(TCa, TCa_apex[i], 'TCa_ZERO_exclude_apex.ipacti',
                                 'ForwardSolveActivation/TCa_'+str(es)+'.ipacti')
        copy('ForwardSolveActivation/TCa_'+str(es)+'.ipacti', 'TCa_current.ipacti')
        copy('OptimisedActivation/TCa_'+str(ed)+'.ipacti', 'TCa_previous.ipacti')

        os.system('cm SolveInitialSystoleMF.com')
        print 'LOG: Forward solution has finished and is converged.\n'
        time.sleep(5)
        active_save_solutions(str(es), 'ForwardSolveExfile', 'ForwardSolveError', 'ForwardSolveCavityVolume',
                              'ForwardSolveStressStrain', 0)

        # b) Set sensible upper and lower bounds for TCa parameter according to pressure change.
        # Get bounds for TCa
        [lb, ub] = active_bounds_get(pressure, int(ed), int(es))

        # c) Optimise TCa parameter for current frame to minimise fitting error at current frame.
        if TCa_apex[i] > 0:
            apex_acti = 1
        else:
            apex_acti = 0
        TCa_Opt = fmin_l_bfgs_b(optimise_active_obj_function, [TCa], args=[study_id, str(es), apex_acti], approx_grad=1,
                                bounds=[(lb, ub)], factr=1e12, epsilon=0.5)

        print 'TCa_Opt: ' + str(TCa_Opt)
        copy('OptimisedExfile/LVContraction_'+str(es)+'.exnode', 'OptimisedExfile/LVContraction_'+str(es)+'_apex_'
             + str(i) + '.exnode')
        # d) Evaluate final optimised solution and fitting errors.
        TCa = float(TCa_Opt[0][0])
        mse[i] = optimise_active_obj_function(TCa, study_id, str(es), apex_acti)

        if i == 0:
            print 'LOG_STATUS_CURRENT: Optimised Passive Parameter with apical activation: ' + str(TCa) + ' for frame '+str(es)
            print 'LOG_STATUS_CURRENT: With total MSE of with apical activation: ' + str(mse[i])
        else:
            print 'LOG_STATUS_CURRENT: Optimised Passive Parameter without apical activation: ' + str(TCa) + ' for frame '+str(es)
            print 'LOG_STATUS_CURRENT: With total MSE of without apical activation: ' + str(mse[i])

    print 'LOG: The MSE values for the two apical contraction methods are: \n'
    print 'LOG: With apical contraction: ', str(mse[0])
    print 'LOG: Without apical contraction: ', str(mse[1])

    if mse[0] > mse[1]:
        apex_acti = 0
        print 'LOG_STATUS: Method excluding apical activation gives better fit at ES.\n'
        time.sleep(5)
    else:
        apex_acti = 1
        print 'LOG_STATUS: Method using apical activation gives better fit at ES.\n'
        time.sleep(5)
    print 'LOG_STATUS: Completed optimise_determine_apical_activation.\n'
    return apex_acti


def optimised_passive_displacement_mse(study_id, study_frame):
    os.chdir(os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/PassiveMechanics')
    C1_opt = 0
    C1_opt = material_get(os.environ['RESULTS']+study_id+'/LV_CubicGuc_Optimised.ipmate')
    copy(os.environ['RESULTS'] + study_id + '/LV_CubicGuc_Optimised.ipmate', 'LV_CubicGuc_Opt.ipmate')

    # Evaluate displacement mse using high stiffness value.
    idx = passive_loop_index(study_frame)
    print 'LOG_STATUS: Evaluating displacement at each frame by projecting surface data onto high stiffness models.'
    time.sleep(5)
    displacement_mse = optimise_passive_obj_function(15, study_id, study_frame, toggle=2)
    os.system('cp -r OptimisedExfile/ DisplacementExfile/')
    with open('DisplacementMSE_high_stiffness.txt', 'w') as fw:
        fw.write('Frame\tMSE (mm^2)\n')
        for i in range(1, len(idx)):
            fw.write(str(idx[i])+'\t'+str(displacement_mse[i-1])+'\n')
        mse_tot = 0
        for i in range(1, len(idx)):
            mse_tot += displacement_mse[i-1]
        fw.write('Total mse\t'+str(mse_tot)+'\n')
    # Restore files to optimal C1.
    print 'LOG_STATUS: Restore files to optimal C1. '
    mse = optimise_passive_obj_function(C1_opt, study_id, study_frame, toggle=1)

    print 'LOG_STATUS: Finished evaluating displacement using surface data projections'


def optimise_passive_parameter_sweep(study_id, study_frame):
    print 'LOG: Doing parameter sweep for C1 parameter of study ' + str(study_id) + '...'
    os.chdir(os.environ['STUDIES'] + study_id + '/LVMechanics' + study_id + '/PassiveMechanics')
    f = open('ParameterSweep.txt', 'a')
    idx = passive_loop_index(study_frame)
    header = ["C1(kPa)"]
    for i in range(1, len(idx)):
        header = [str(header[0]) + "\tFrame " + str(idx[i])]
    header = [str(header[
                      0]) + "\tTotal mse (mm^2)\tED Equatorial mean strain\Strain std\ED Equatorial total stress\tStress std\n"]

    copy(os.environ['RESULTS'] + study_id + '/LV_CubicGuc_Optimised.ipmate', 'LV_CubicGuc_Opt.ipmate')

    copy('LV_CubicGuc_Opt.ipmate', 'LV_CubicGuc.ipmate')
    C1_opt = material_get('LV_CubicGuc.ipmate')
    f.write(str(header[0]))
    # C1_sweep = [C1_opt, 0.8, 0.9, 1.0, C1_opt]
    start = numpy.floor(C1_opt)
    if start < 1:
        start = 1
    # start = numpy.ceil(C1_opt)
    # C1_sweep = numpy.linspace(start, 10, (10-start+1))
    C1_sweep = numpy.linspace(2, 10, 9)
    # C1_sweep = numpy.insert(C1_sweep, 0, [start - 0.2])
    C1_sweep = numpy.append(C1_sweep, C1_opt)
    # C1_sweep = [start-0.5, start-0.2, C1_sweep, C1_opt]
    # C1_sweep = [C1_opt-0.1] +  C1_sweep +  [C1_opt]
    for i in range(0, len(C1_sweep)):
        # 1. Update current C1 estimate in ipmate file.
        copy('LV_CubicGuc.ipmate', 'LV_CubicGuc_previous.ipmate')
        material_create_ipmate(C1_sweep[i], 'LV_CubicGuc_TEMPLATE.ipmate', 'LV_CubicGuc.ipmate')
        # """
        for j in range(1, len(idx)):
            print C1_sweep
            print 'LOG: Evaluating C1 as ', str(C1_sweep[i]), ' for frame ', str(idx[j])
            time.sleep(5)
            # Solve warm start using new C1 guess - using non-registration of surface data.
            passive_warm_solve(study_id, str(idx[j]))
            quit()
        # """
        # 2. Evaluate mean equatorial stress and strain.
        [strain, strain_std, stress, stress_std] = _evaluate_equatorial_stress(study_id, idx)
        # 3. Evaluate MSE of fitting and sum for all frames from DS+1 to ED.
        mse_ = []
        for j in range(1, len(idx)):
            # Vector error value
            info = open('OptimisedError/EpiError_' + str(idx[j]) + '.opdata').read()
            array_ = re.split(r'[\n\\]', info)
            temp = array_[3].split()
            num_epi = float(temp[len(temp) - 1])
            temp = array_[7].split()
            epi_rmse = float(temp[len(temp) - 1])

            info = open('OptimisedError/EndoError_' + str(idx[j]) + '.opdata').read()
            array_ = re.split(r'[\n\\]', info)
            temp = array_[3].split()
            num_endo = float(temp[len(temp) - 1])
            temp = array_[7].split()
            endo_rmse = float(temp[len(temp) - 1])

            mse_.append((epi_rmse ** 2 * num_epi + endo_rmse ** 2 * num_endo) / (num_endo + num_epi))
            print '\033[0;30;45m LOG: Current MSE for frame ' + str(idx[j]) + ' = ' + str(mse_[j - 1]) + '\033[0m\n'
        mse_tot = 0
        line = [str(C1_sweep[i])]
        for j in range(1, len(idx)):
            mse_tot = mse_tot + mse_[j - 1]
            line = [str(line[0]) + "\t" + str(mse_[j - 1])]
        line = [str(line[0]) + "\t" + str(mse_tot) + "\t" + str(strain[-1]) + "\t" + str(strain_std[-1]) + "\t" +
                str(stress[-1]) + "\t" + str(stress_std[-1]) + "\n"]
        f.write(str(line[0]))
        print 'LOG: Current total MSE for diastole = ' + str(mse_tot) + '\n'

    print 'LOG: Finished parameter sweep on C1'
    f.close()


def optimise_passive_parameter_sweep_ed_only(study_id, study_frame):
    print 'LOG: Doing parameter sweep for C1 parameter of study ' + str(study_id) + '...'
    os.chdir(os.environ['STUDIES'] + study_id + '/LVMechanics' + study_id + '/PassiveMechanics')
    if not os.path.exists('ParameterSweep_ED'):
        os.mkdir('ParameterSweep_ED')
    f = open('ParameterSweep_ED.txt', 'a')
    fstress = open('EquatorialStress_ED.txt', 'a')
    fstress.write('ED Equatorial mean strain\tStrain std\tED Equatorial total stress\tStress std\tGlobal circ strain\tstrain std\tGlobal long strain\tstrain std\tGlobal radial strain\tstrain std\n')
    idx = passive_loop_index(study_frame)
    header = ["C1(kPa)"]
    for i in range(0, 16):
        header = [str(header[0]) + "\tEpi " + str(i+1)]
    for i in range(0, 16):
        header = [str(header[0]) + "\tEndo " + str(i + 1)]
    header = [str(header[0]) + '\tEpi MSE\tEndo MSE\tTotal MSE\tSigned Epi MSE\tSigned Endo MSE\tSigned total MSE']
    header = [str(header[0])+'\n']
    copy(os.environ['RESULTS'] + study_id + '/LV_CubicGuc_Optimised.ipmate', 'LV_CubicGuc_Opt.ipmate')

    copy('LV_CubicGuc_Opt.ipmate', 'LV_CubicGuc.ipmate')
    C1_opt = material_get('LV_CubicGuc.ipmate')
    f.write(str(header[0]))
    start = numpy.floor(C1_opt)
    if start < 1:
        start = 1
    #C1_sweep = numpy.linspace(2 10, 9)
    if (C1_opt >= 3) & (C1_opt < 7):
        C1_sweep = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
    elif (C1_opt >= 2) & (C1_opt < 3):
        C1_sweep = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
    elif (C1_opt >= 1) & (C1_opt < 2):
        C1_sweep = [0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
    elif C1_opt < 1:
        C1_sweep = [0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
    elif C1_opt >= 7:
        C1_sweep = [3, 4, 5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
    else:
        print 'Log fail for parameter sweep.'
        quit()
    if (study_id == 'STF_18') | (study_id == 'STF_16') | (study_id == 'STF_20'):
        #C1_sweep = [0.4, 0.6, 0.8]
        C1_sweep = [0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
    elif study_id == 'MR_087087':
        C1_sweep = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
    elif study_id == 'MR_160160':
        C1_sweep = [1.0, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
    elif study_id == 'MR_124124':
        C1_sweep = [1.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
    C1_sweep = numpy.append(C1_sweep, C1_opt)
    regional_error_epi = numpy.zeros([len(C1_sweep), 16])
    regional_error_endo = numpy.zeros([len(C1_sweep), 16])
    signed_regional_error_endo = numpy.zeros([len(C1_sweep), 16])
    signed_regional_error_epi =  numpy.zeros([len(C1_sweep), 16])
    data_points_excluded_epi = numpy.loadtxt('ForwardSolveError/DataPointsExcludedEpi.txt', float)
    data_points_excluded_endo = numpy.loadtxt('ForwardSolveError/DataPointsExcludedEndo.txt', float)
    for i in range(0, len(C1_sweep)):
        # 1. Update current C1 estimate in ipmate file.
        copy('LV_CubicGuc.ipmate', 'LV_CubicGuc_previous.ipmate')
        material_create_ipmate(C1_sweep[i], 'LV_CubicGuc_TEMPLATE.ipmate', 'LV_CubicGuc.ipmate')
        print 'LOG: Evaluating C1 as ', str(C1_sweep[i]), ' for ED frame'
        # Solve warm start using new C1 guess - using non-registration of surface data.
        passive_warm_solve(study_id, '1')
        if i < (len(C1_sweep)-1):
            copy('OptimisedError/EndoError_1.exdata', 'ParameterSweep_ED/EndoError_'+str(i)+'.exdata')
            copy('OptimisedError/EpiError_1.exdata', 'ParameterSweep_ED/EpiError_' + str(i) + '.exdata')
            copy('OptimisedExfile/LVInflation_1.exnode', 'ParameterSweep_ED/LVInflationED_'+str(i)+'.exnode')
            copy('OptimisedExfile/LVInflation_1.exelem', 'ParameterSweep_ED/LVInflationED_' + str(i) + '.exelem')
        # 2. Evaluate average error on a regional basis.
        [coord_epi, projection_epi] = _read_exdata_values('OptimisedError/EpiError_1.exdata')
        [coord_endo, projection_endo] = _read_exdata_values('OptimisedError/EndoError_1.exdata')
        [elem_epi, xi_epi] = _read_xi_values('OptimisedError/ProjectionXiEpi_1.ipxi')
        [elem_endo, xi_endo] = _read_xi_values('OptimisedError/ProjectionXiEndo_1.ipxi')
        [elem_epi_contain, xi_epi_contain] = _read_xi_values('OptimisedError/ContainXiEpi_1.ipxi')
        [elem_endo_contain, xi_endo_contain] = _read_xi_values('OptimisedError/ContainXiEndo_1.ipxi')
        TOL = 1e-2
        elem_counter_epi = numpy.zeros(16)
        elem_counter_endo = numpy.zeros(16)
        for j in range(0, len(coord_endo)):
            if not j in data_points_excluded_epi:
                if elem_epi_contain[j] == 0:
                    signed_error = -numpy.linalg.norm(projection_epi[j])**2 # If data point is outside of the model, the error is negative.
                else:
                    signed_error = numpy.linalg.norm(projection_epi[j])**2 # If data point is inside the model, the error is positive.
                regional_error_epi[i][elem_epi[j]-1] = regional_error_epi[i][elem_epi[j]-1] + numpy.linalg.norm(projection_epi[j])**2
                signed_regional_error_epi[i][elem_epi[j]-1] = signed_regional_error_epi[i][elem_epi[j]-1] + signed_error
                elem_counter_epi[elem_epi[j]-1] = elem_counter_epi[elem_epi[j]-1] + 1
            if not j in data_points_excluded_endo:
                if elem_endo_contain[j] == 0:
                    signed_error = numpy.linalg.norm(projection_endo[j])**2 # If data point is outside of the model, the error is positive.
                else:
                    signed_error = -numpy.linalg.norm(projection_endo[j])**2 # If data point is inside the model, the error is negative.
                regional_error_endo[i][elem_endo[j] - 1] = regional_error_endo[i][elem_endo[j] - 1] + numpy.linalg.norm(projection_endo[j])**2
                signed_regional_error_endo[i][elem_endo[j] - 1] = signed_regional_error_endo[i][elem_endo[j] - 1] + signed_error
                elem_counter_endo[elem_endo[j] - 1] = elem_counter_endo[elem_endo[j] - 1] + 1
        # Average the error per region
        for j in range(0, 16):
            regional_error_epi[i][j] = regional_error_epi[i][j] / elem_counter_epi[j]
            regional_error_endo[i][j] = regional_error_endo[i][j] / elem_counter_endo[j]
            signed_regional_error_epi[i][j] = signed_regional_error_epi[i][j] / elem_counter_epi[j]
            signed_regional_error_endo[i][j] = signed_regional_error_endo[i][j] / elem_counter_endo[j]
        line = [str(C1_sweep[i])]
        sum_epi = 0
        sum_endo = 0
        signed_sum_epi = 0
        signed_sum_endo = 0
        for j in range(0, 16):
            line = [str(line[0]) + "\t" + str(regional_error_epi[i][j])]
            sum_epi = sum_epi + regional_error_epi[i][j]
            signed_sum_epi = signed_sum_epi + signed_regional_error_epi[j][j]
        for j in range(0, 16):
            line = [str(line[0]) + "\t" + str(regional_error_endo[i][j])]
            sum_endo = sum_endo + regional_error_endo[i][j]
            signed_sum_endo = signed_sum_endo + signed_regional_error_endo[j][j]
        # Write total MSE
        line = [str(line[0]) + "\t" + str(sum_epi/16.0) + "\t" +
        str(sum_endo/16.0) + "\t" + str((sum_epi+sum_endo)/32.0) +
        "\t" + str(signed_sum_epi/16.0)+ "\t" + str(signed_sum_endo/16.0)+
         "\t" + str((signed_sum_epi+signed_sum_endo)/32.0)]
        line = [str(line[0])+"\n"]
        f.write(str(line[0]))
        # Evaluate equatorial stress
        ed_idx = '1';
        output = _evaluate_equatorial_stress(study_id, ed_idx)
        line = [str(C1_sweep[i])+"\t"+str(output[0]) + "\t" +
                str(output[1]) + "\t" + str(output[2]) + "\t" +
                str(output[3]) + "\t" + str(output[4])+ "\t" +
                str(output[5]) + "\t" + str(output[6])+ "\t" +
                str(output[7]) + "\t" + str(output[8]) + "\t" +
                str(output[9]) + "\n"]
        fstress.write(line[0])


def _evaluate_equatorial_stress(study_id, idx):
    os.chdir(os.environ['STUDIES'] + study_id + '/LVMechanics' + study_id + '/PassiveMechanics')
    #for i in range(1, len(idx)):
    filename = 'OptimisedStressStrain/Strain_' + str(idx) + '.exdata'
    [elem, xi, strain] = _read_exdata_stress_values(filename)
    filename = 'OptimisedStressStrain/Strain_wall_' + str(idx) + '.exdata'
    [elem_wall, xi_wall, strain_wall] = _read_exdata_stress_values(filename)
    filename = 'OptimisedStressStrain/TotalStress_' + str(idx) + '.exdata'
    [elem, xi, stress] = _read_exdata_stress_values(filename)
    elem_strain_sum = []
    elem_stress_sum = []
    radial_strain_sum = []
    circ_strain_sum = []
    long_strain_sum = []
    radial_stress_sum = []
    circ_stress_sum = []
    long_stress_sum = []
    for i in range(0, len(elem)):
        circ_strain_sum.append(strain_wall[i][0])
        long_strain_sum.append(strain_wall[i][1])
        radial_strain_sum.append(strain_wall[i][2])
        if (elem[i] < 12) & (elem[i] > 7):
            elem_strain_sum.append(strain[i][0])
            elem_stress_sum.append(stress[i][0])
    elem_strain_sum = numpy.array(elem_strain_sum)
    elem_stress_sum = numpy.array(elem_stress_sum)
    strain_mean = numpy.mean(elem_strain_sum)
    strain_std = numpy.std(elem_strain_sum)
    stress_mean = numpy.mean(elem_stress_sum)
    stress_std = numpy.std(elem_stress_sum)
    circ_strain_mean = numpy.mean(circ_strain_sum)
    circ_strain_std = numpy.std(circ_strain_sum)
    long_strain_mean = numpy.mean(long_strain_sum)
    long_strain_std = numpy.std(long_strain_sum)
    radial_strain_mean = numpy.mean(radial_strain_sum)
    radial_strain_std = numpy.std(radial_strain_sum)

    output = [strain_mean, strain_std, stress_mean, stress_std, circ_strain_mean, circ_strain_std, long_strain_mean, long_strain_std, radial_strain_mean, radial_strain_std]
    return output


def optimise_passive_parameter_sweep_sensible_projection(study_id, study_frame):
    print 'LOG: Doing parameter sweep for C1 parameter of study ' + str(study_id) + '...'
    os.chdir(os.environ['STUDIES'] + study_id + '/LVMechanics' + study_id + '/PassiveMechanics')
    f = open('ParameterSweep.txt', 'a')
    idx = passive_loop_index(study_frame)
    header = ["C1(kPa)"]
    for i in range(1, len(idx)):
        header = [str(header[0]) + "\tFrame " + str(idx[i])]
    header = [str(header[
                      0]) + "\tTotal mse (mm^2)\tED Equatorial mean strain\Strain std\ED Equatorial total stress\tStress std\n"]

    copy(os.environ['RESULTS'] + study_id + '/LV_CubicGuc_Optimised.ipmate', 'LV_CubicGuc_Opt.ipmate')

    copy('LV_CubicGuc_Opt.ipmate', 'LV_CubicGuc.ipmate')
    C1_opt = material_get('LV_CubicGuc.ipmate')
    f.write(str(header[0]))
    # C1_sweep = [C1_opt, 0.8, 0.9, 1.0, C1_opt]
    start = numpy.floor(C1_opt)
    if start < 1:
        start = 1
    # start = numpy.ceil(C1_opt)
    # C1_sweep = numpy.linspace(start, 10, (10-start+1))
    C1_sweep = numpy.linspace(2, 10, 9)
    # C1_sweep = numpy.insert(C1_sweep, 0, [start - 0.2])
    C1_sweep = numpy.append(C1_sweep, C1_opt)
    # C1_sweep = [start-0.5, start-0.2, C1_sweep, C1_opt]
    # C1_sweep = [C1_opt-0.1] +  C1_sweep +  [C1_opt]
    for i in range(0, len(C1_sweep)):
        # 1. Update current C1 estimate in ipmate file.
        copy('LV_CubicGuc.ipmate', 'LV_CubicGuc_previous.ipmate')
        material_create_ipmate(C1_sweep[i], 'LV_CubicGuc_TEMPLATE.ipmate', 'LV_CubicGuc.ipmate')
        for j in range(1, len(idx)):
            print C1_sweep
            print 'LOG: Evaluating C1 as ', str(C1_sweep[i]), ' for frame ', str(idx[j])
            time.sleep(5)
            # Solve warm start using new C1 guess - using non-registration of surface data.
            passive_warm_solve(study_id, str(idx[j]))
        # 2. Evaluate mean equatorial stress and strain.
        [strain, strain_std, stress, stress_std] = _evaluate_equatorial_stress(study_id, idx)
        # 3. Evaluate MSE of fitting and sum for all frames from DS+1 to ED.
        error_magnitude = optimise_passive_obj_evaluate_sensible_projection(idx, 2)
        total_error_magnitude = 0
        line = [str(C1_sweep[i])]
        for j in range(1, len(idx)):
            total_error_magnitude = total_error_magnitude + error_magnitude[j - 1]
            line = [str(line[0]) + "\t" + str(error_magnitude[j - 1])]
        line = [str(line[0]) + "\t" + str(total_error_magnitude) + "\t" + str(strain[-1]) + "\t" + str(
            strain_std[-1]) + "\t" +
                str(stress[-1]) + "\t" + str(stress_std[-1]) + "\n"]
        f.write(str(line[0]))

    print 'LOG: Finished parameter sweep on C1'
    f.close()
