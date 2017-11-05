from passive_mechanics import *
from optimise import optimise_passive_obj_function
import os, sys

def setup_get_frames(idx):
    # This function gets the specified study ID and important frame numbers.
    filename = os.environ['PARAM_ESTIMATION'] + '/NYStFranFrameNumber_UsedForAnalysis.txt'
    f = open(filename, 'r')

    study_ids = []
    study_frames = []
    c1_inits = []
    study_info = f.readline()
    while len(study_info) != 0:  # Reaching the end of the file
        study_ids.append(study_info.split()[0])
        study_frames.append(study_info.split()[1:5])
        c1_inits.append(study_info.split()[5])
        study_info = f.readline()
        num_studies = len(study_ids)

    study_id = study_ids[idx]
    study_frame = study_frames[idx]
    c1_init = float(c1_inits[idx])

    return study_id, study_frame, c1_init


def hessian_convergence(study_name, study_frame):
    # Main commands begin here.
    #os.system('rm *.*~')  # Remove all temporary files in Main folder.

    # Top level control of the entire simulation:
    # Command line inputs:
    epsilons = [0.01, 0.05, 0.08, 0.1]
    #epsilons = [0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2]
    #epsilons = [0.6, 0.4, 0.2, 0.1, 0.08, 0.06, 0.04, 0.02, 0.01, 0.008, 0.006, 0.004, 0.002, 0.001]
    with open('HessianConvergence_' + study_name + '.txt', 'a') as f:
        for i in range(0, len(epsilons)):
            eps = epsilons[i]
            d2Error = passive_identifiability_evaluate(study_name, study_frame, eps)
            f.write(str(eps) + '\t' + str(d2Error) + '\n')


def hessian_test_epsilon(c1_opt, opt_val, study_name, study_frame):
    study_dir = os.environ['STUDIES'] + study_name
    os.chdir(study_dir)
    idx = passive_loop_index(study_frame)
    mse_opt = opt_val
    #mse_opt = optimise_passive_obj_function(c1_opt, study_name, study_frame, 1)
    with open('HessianConvergence_' + study_name + '.txt', 'a+') as f:
        f.write('Epsilon\tChange in MSE produced\n')
        for i in range(0, 100):
            eps = 0.05 * (i + 1)
            mse_plus = optimise_passive_obj_function((c1_opt + eps), study_name, study_frame, 1)
            mse_minus = optimise_passive_obj_function((c1_opt - eps), study_name, study_frame, 1)
            d2Error = 1 / (eps * eps) * (mse_plus - 2 * mse_opt + mse_minus)
            f.write(str(eps) + '\t' + str(abs(mse_plus - mse_opt))+ '\t' + str(abs(mse_minus - mse_opt))+'\t' + str(d2Error)+'\n')
            if min([abs(mse_plus - mse_opt), abs(mse_minus - mse_opt)]) > 0.01*len(idx):
                f.write('Smallest epsilon allowed\tMSE plus\tMSE optimal\tMSE minus\tSecond derivative calculation\n')
                #d2Error = passive_identifiability_evaluate(study_name, study_frame, eps)
                d2Error = 1 / (eps * eps) * (mse_plus - 2 * mse_opt + mse_minus)
                f.write(str(eps) + '\t'  + str(mse_plus) + '\t' + str(mse_opt) + '\t' + str(
                    mse_minus)+ '\t' + str(d2Error) + '\n')
                print 'Smallest epsilon allowed | Second derivative:\n'
                print str(eps) + '\t' + str(d2Error) + '\n'
                print 'MSE plus | Optimal MSE | MSE minus\n'
                print str(mse_plus) + '\t' + str(mse_opt) + '\t' + str(mse_minus) + '\n'
                break
    print 'Finished 2nd derivative calculation.'
    return eps, d2Error, len(idx)

if __name__ == '__main__':
    idx = int(sys.argv[1])  # Study number
    [study_name, study_frame, c1_init] = setup_get_frames(idx)
    hessian_convergence(study_name, study_frame)
