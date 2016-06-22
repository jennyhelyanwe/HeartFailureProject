__author__ = 'zwan145'

import os
from material_setup import *
from mech_setup import *
from bc_setup import *
import fileinput


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
    copy('LV_CubicPreEpiBase_TEMPLATE.ipinit', 'ForwardSolveSolution/CurrentInflated_'+str(ds)+'.ipinit')

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
        bc_pressure_set(p_increm, 'ForwardSolveSolution/CurrentInflated_'+str(idx[i])+'.ipinit',
                        'LV_CubicPreEpiBase.ipinit')

        # Update displacement BC
        data_cur_epi = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Epi_'+str(idx[i])+'.ipdata'
        data_next_epi = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Epi_'+str(idx[i+1])+'.ipdata'

        if i == 0:
            # Get basal nodal coordinates at DS
            filename = 'DSWall.ipnode'
        else:
            # Get basal nodal coordinates at previous frame.
            filename = 'ForwardSolveExfile/LVInflation_'+str(idx[i])+'.exnode'
        data_cur_endo = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Endo_'+str(idx[i])+'.ipdata'
        data_next_endo = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Endo_'+str(idx[i+1])+'.ipdata'
        bc_displacement_set(node_idx, data_cur_epi, data_next_epi, data_cur_endo, data_next_endo,
                            'LV_CubicPreEpiBase.ipinit', filename)

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
    os.system('cm SolveWarmStartPassive_nonReg.com')


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
    copy('output/LVInflation_gauss_ER.exdata', stress_folder+'/ER_'+frame_num+'.exdata')
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
        ## 5a. Save registered surface data for visualisation in CMGUI.
        copy('output_debug/ED_Endo_nonreg.exdata', 'OptimisedExfile/Surface_endo_'+frame_num+'.exdata')
        copy('output_debug/ED_Epi_nonreg.exdata', 'OptimisedExfile/Surface_epi_'+frame_num+'.exdata')
        ## 6. Save current cavity volumes.
        copy('output_cavity_volume/LVCavityUpdate.opelem', cavity_folder+'/LVCavity_'+frame_num+'.opelem')
    else:
        ## 5. Save warm-start model prediction/solutions.
        copy('output/LVInflation.ipinit', 'ForwardSolveSolution/CurrentInflated_'+frame_num+'.ipinit')
        copy('LV_CubicGuc.ipmate', 'LV_CubicGuc_ForwardSolve.ipmate')
        copy('output_debug/ED_Endo_nonreg.exdata', 'ForwardSolveExfile/Surface_endo_'+frame_num+'.exdata')
        copy('output_debug/ED_Epi_nonreg.exdata', 'ForwardSolveExfile/Surface_epi_'+frame_num+'.exdata')
        ## 6. Save current cavity volumes.
        copy('output_cavity_volume/LVCavityCurrent.opelem', cavity_folder+'/LVCavity_'+frame_num+'.opelem')
#
#=======================================================================================================================
#

def passive_visualisation_customise(study_id, study_frame):
    ds, ed, es, tot = tuple(study_frame)
    os.chdir(os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/PassiveMechanics')
    # Modify visualisation file to be specific to study frames of current study.
    with open('ViewDisplacementMSE.com', 'r') as f:
        filedata = f.read()
    tempdata = filedata.replace("$DS = 23", "$DS = "+ds)
    newdata = tempdata.replace("$tot = 30", "$tot = "+tot)
    with open('ViewDisplacementMSE.com', 'w') as f:
        f.write(newdata)

    # Modify visualisation script for forward solve solutions.
    with open('ViewForwardSolve.com', 'r') as f:
        filedata = f.read()
    tempdata = filedata.replace("$DS = 23", "$DS = "+ds)
    newdata = tempdata.replace("$tot = 30", "$tot = "+tot)
    with open('ViewForwardSolve.com', 'w') as f:
        f.write(newdata)

    # Modify visualisation of optimised solution.
    with open('ViewOptimisedModels.com', 'r') as f:
        filedata = f.read()
    tempdata = filedata.replace("$DS = 23", "$DS = "+ds)
    newdata = tempdata.replace("$tot = 30", "$tot = "+tot)
    with open('ViewOptimisedModels.com', 'w') as f:
        f.write(newdata)


def passive_displacement_mse(study_id, study_frame):
    # This function projects each surface data frame in rapid filling diastole onto the diastasis geometry to quantify
    # the projected geometric displacement in mean squared error format. Then writes this information into a text file.
    ds, ed, es, tot = tuple(study_frame)

    os.chdir(os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/PassiveMechanics')

    if not os.path.exists('DisplacementMSE'):
        os.mkdir('DisplacementMSE')

    # Loop through all frames in passive portion of simulation
    idx = passive_loop_index(study_frame)
    for i in range(0, len(idx)):
        fname_epi = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Epi_'+str(idx[i])
        fname_endo = os.environ['GEOM_DATA']+study_id+'/Passive/'+study_id+'_Surface_Points_Endo_'+str(idx[i])
        fname_epi_error = 'DisplacementMSE/EpiDisp_'+str(idx[i])
        fname_endo_error = 'DisplacementMSE/EndoDisp_'+str(idx[i])
        fname_epi_data = 'DisplacementMSE/EpiData_'+str(idx[i])
        fname_endo_data = 'DisplacementMSE/EndoData_'+str(idx[i])
        # Modify template command file with directory.
        with open('DisplacementMSE_TEMPLATE.com','r') as f:
            with open('DisplacementMSE.com','w') as f_w:
                temp = f.readline()
                while temp != '':
                    if temp.find('TEMPLATE_EPI') >= 0:
                        f_w.write('fem def data;r;'+fname_epi+'\n')
                    elif temp.find('TEMPLATE_DATA_EPI') >= 0:
                        f_w.write('fem exp data;'+fname_epi_data+' as Epi_data\n')
                    elif temp.find('TEMPLATE_ERROR_EPI') >= 0:
                        f_w.write('fem list data;'+fname_epi_error+' error\n')
                    elif temp.find('TEMPLATE_EX_ERROR_EPI') >= 0:
                        f_w.write('fem exp data;'+fname_epi_error+' as Epi_error error')
                    elif temp.find('TEMPLATE_ENDO') >= 0:
                        f_w.write('fem def data;r;'+fname_endo+'\n')
                    elif temp.find('TEMPLATE_DATA_ENDO') >= 0:
                        f_w.write('fem exp data;'+fname_endo_data+' as Endo_data\n')
                    elif temp.find('TEMPLATE_ERROR_ENDO') >= 0:
                        f_w.write('fem list data;'+fname_endo_error+' error\n')
                    elif temp.find('TEMPLATE_EX_ERROR_ENDO') >= 0:
                        f_w.write('fem exp data;'+fname_endo_error+' as Endo_error error')
                    else:
                        f_w.write(temp)
                    temp = f.readline()

        os.system('cm DisplacementMSE.com')

    with open('DisplacementMSE/Summary.txt','w') as f:
        f.write('Frame\tEpi disp.\tEndo disp.\tAverage disp.\n')
        sum_mse = 0
        for i in range(0, len(idx)):
            info = open('DisplacementMSE/EpiDisp_'+str(idx[i])+'.opdata').read()
            array = re.split(r'[\n\\]', info)
            temp = array[3].split()
            num_epi = float(temp[len(temp)-1])
            temp = array[7].split()
            epi_rmse = float(temp[len(temp)-1])
            epi_mse = epi_rmse**2

            info = open('DisplacementMSE/EndoDisp_'+str(idx[i])+'.opdata').read()
            array = re.split(r'[\n\\]', info)
            temp = array[3].split()
            num_endo = float(temp[len(temp)-1])
            temp = array[7].split()
            endo_rmse = float(temp[len(temp)-1])
            endo_mse = endo_rmse**2

            ave_mse = (epi_rmse**2*num_epi + endo_rmse**2*num_endo)/(num_endo+num_epi)
            sum_mse = sum_mse + ave_mse
            f.write(str(idx[i])+'\t'+str(endo_mse)+'\t'+str(epi_mse)+'\t'+str(ave_mse)+'\n')
        f.write('Total average MSE over all frames:\n')
        f.write(str(sum_mse))


def passive_identifiability_evaluate(study_id, study_frame):
    # This function evaluates the sensitivity to the TCa parameter at each frame by perturbing the TCa parameter
    # and evaluating the 1st and 2nd derivatives for the resulting fitting error.
    print 'LOG: Evaluating C1 identifiability (2nd derivative of objective function)...'
    os.chdir(os.environ['STUDIES'] + study_id + '/LVMechanics' + study_id + '/PassiveMechanics')
    f = open('StiffnessIdentifiability.txt', 'w')
    f.write('C1\t1st_d\t2nd_d\tmse(x+2h)\tmse(x+h)\tmse(x)\tmse(x-h)\tmse(x-2h)\n')
    copy(os.environ['RESULTS'] + study_id + '/LV_CubicGuc_Optimised.ipmate', 'LV_CubicGuc_Opt.ipmate')

    copy('LV_CubicGuc_Opt.ipmate', 'LV_CubicGuc.ipmate')
    C1_opt = material_get('LV_CubicGuc.ipmate')
    C1 = C1_opt
    print C1
    eps = 0.1
    h = numpy.array([eps*2, eps, -eps, -eps*2, 0.0])  # Evaluate opt C1 last so that the folders are restored to opt state.
    assert isinstance(C1, object)
    test = h + C1
    print test
    C1 = [test[0], test[1], test[4], test[2], test[3]]  # Re-order for derivative calculations.
    print 'C1 values evaluated are: \n', str(C1)
    time.sleep(10)
    if min(test) > 0:
        mse = [0, 0, 0, 0, 0]  # Space holder mse values.
        for j in range(0, len(h)):
            # 1. Update current C1 estimate in ipmate file.
            copy('LV_CubicGuc.ipmate', 'LV_CubicGuc_previous.ipmate')
            material_create_ipmate(test[j], 'LV_CubicGuc_TEMPLATE.ipmate', 'LV_CubicGuc.ipmate')

            # 2. Re-run simulation using current C1 estimate and get new model predictions for
            # each frame from DS+1 to ED.
            idx = passive_loop_index(study_frame)
            for i in range(1, len(idx)):
                # Solve warm start using new C1 guess - using non-registration of surface data.
                passive_warm_solve(study_id, str(idx[i]))

            # 3. Evaluate MSE of fitting and sum for all frames from DS+1 to ED.
            mse_ = []
            for i in range(1, len(idx)):
                # Vector error value
                info = open('OptimisedError/EpiError_' + str(idx[i]) + '.opdata').read()
                array_ = re.split(r'[\n\\]', info)
                temp = array_[3].split()
                num_epi = float(temp[len(temp) - 1])
                temp = array_[7].split()
                epi_rmse = float(temp[len(temp) - 1])

                info = open('OptimisedError/EndoError_' + str(idx[i]) + '.opdata').read()
                array_ = re.split(r'[\n\\]', info)
                temp = array_[3].split()
                num_endo = float(temp[len(temp) - 1])
                temp = array_[7].split()
                endo_rmse = float(temp[len(temp) - 1])

                mse_.append((epi_rmse ** 2 * num_epi + endo_rmse ** 2 * num_endo) / (num_endo + num_epi))
                print '\033[0;30;45m LOG: Current MSE for frame ' + str(idx[i]) + ' = ' + str(mse_[i - 1]) + '\033[0m\n'
            mse_tot = 0
            for i in range(1, len(idx)):
                mse_tot = mse_tot + mse_[i - 1]
            print 'LOG: Current total MSE for diastole = ' + str(mse_tot) + '\n'
            mse[j] = mse_tot
        mse = [mse[0], mse[1], mse[4], mse[2], mse[3]]  # Re-order for derivative calculations.
        C1 = [test[0], test[1], test[4], test[2], test[3]]  # Re-order for derivative calculations.
        print 'C1 values evaluated are: ', str(C1)
        # Get 1st derivative value
        dError = 0.5 * (mse[1] - mse[3])/h[1]
        print 'First derivative at optimised C1: ', str(dError)
        d2Error = 1 / (12 * h[1] * h[1]) * (-mse[0] + 16 * mse[1] - 30 * mse[2] + 16 * mse[3] - mse[4])
        print 'Second derivative at optimised C1: ', str(d2Error)
        f.write(
            str(C1_opt) + '\t' + str(dError) + '\t' + str(d2Error) + '\t' + str(mse[0]) + '\t' + str(mse[1]) + '\t' +
            str(mse[2]) + '\t' + str(mse[3]) + '\t' + str(mse[4]) + '\n')
        print mse
    else:
        print 'LOG: C1 value too low to evaluate derivatives for. '
    copy('LV_CubicGuc_Opt.ipmate', 'LV_CubicGuc.ipmate')
    f.close()
    print 'LOG: Finished sensitivity analysis.'


def passive_parameter_sweep(study_id, study_frame):
    print 'LOG: Doing parameter sweep for C1 parameter...'
    os.chdir(os.environ['STUDIES'] + study_id + '/LVMechanics' + study_id + '/PassiveMechanics')
    f = open('ParameterSweep.txt', 'w+')
    idx = passive_loop_index(study_frame)
    header = ["C1(kPa)"]
    for i in range(1, len(idx)):
        header = [str(header[0]) + "\tFrame " + str(idx[i])]
    header = [str(header[0]) + "\tTotal mse (mm^2)\n"]


    copy(os.environ['RESULTS'] + study_id + '/LV_CubicGuc_Optimised.ipmate', 'LV_CubicGuc_Opt.ipmate')

    copy('LV_CubicGuc_Opt.ipmate', 'LV_CubicGuc.ipmate')
    C1_opt = material_get('LV_CubicGuc.ipmate')
    f.write(str(header[0]))
    # C1_sweep = [C1_opt, 0.8, 0.9, 1.0, C1_opt]
    start = numpy.floor(C1_opt)
    if start < 1:
        start = 1
    #start = numpy.ceil(C1_opt)
    C1_sweep = numpy.linspace(start, 10, (10-start+1))
    #C1_sweep = numpy.insert(C1_sweep, 0, [start - 0.2])
    C1_sweep = numpy.append(C1_sweep, C1_opt)
    #C1_sweep = [start-0.5, start-0.2, C1_sweep, C1_opt]
    #C1_sweep = [C1_opt-0.1] +  C1_sweep +  [C1_opt]
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
        line = [str(line[0]) + "\t" + str(mse_tot) + "\n"]
        f.write(str(line[0]))
        print 'LOG: Current total MSE for diastole = ' + str(mse_tot) + '\n'

    print 'LOG: Finished parameter sweep on C1'
    f.close()
