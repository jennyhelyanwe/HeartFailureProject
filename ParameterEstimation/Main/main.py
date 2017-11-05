__author__ = 'zwan145'
"""
This main function is designed to automatically implement the following steps
(1) perform passive mechanics simulation and estimate C1
(2) Use optimal C1 to perform active mechanics simulation and estimate Tca
-- Author: Vicky Wang
-- Modified for multi-frame: Jenny Zhinuo Wang
-- Modified: Vicky Wang (26/11/2013)
-- Modified: ZJW (29/01/2015)
-- Modified: ZJW (14/05/2015)
"""

# Import inbuilt subroutines #######
import os
# Import custom python scripts
from setup import *
from geom_setup import *
from results import *
from bc_setup import *
from optimise import *
from passive_mechanics import *
from hessian_convergence import *


def main_program(study_id, study_frame, c1_init, debug_toggle, forward_solve_toggle, a_p_toggle):
    # This function does the following:
    # 1. Program set-up: - set up directories for current study
    #                    - set up log text file if needed.
    # 2. Geometric set-up: - set up reference model and cavity models of LV.
    #                      - transfer reference, cavity models and DS, ED and ES surface data files to study folder.
    # 3. Mechanics simulation set-up: - get pressure array
    #                                 - get data index for prescribing basal displacement boundary conditions.
    #                                 - transfer mechanical simulation template files into study folder.
    # 4. Run estimation: - passive parameter estimation
    #                    - active parameter estimation.
    ####################################################################################################################
    # 1. Program set-up:
    # Set up directories for current study.
    setup_dir(study_id)
    # Set up log text file if needed.

    if debug_toggle == 1:
        if a_p_toggle == 1:
            # Set up log file
            setup_log_file(study_id, 'Output_Passive_sweep.log')
        elif a_p_toggle == 2:
            setup_log_file(study_id, 'Output_Active.log')
        else:
            setup_log_file(study_id, 'Output_P_A.log')

    # 2. Geometric set-up:
    # Set up reference model and cavity model of LV.
    geom_setup_ref_cavity(study_id, study_frame)

    # Initialise mechanics problem with reference model and endpoint surface data.
    geom_setup_data(study_id, study_frame)

    # 3. Mechanics simulation set-up:
    # Get pressure array.
    pressure = bc_pressure_get(study_id, study_frame)

    # Get data indices for applying BC.
    data_idx = bc_displacement_get(study_id)

    # Transfer mechanics template files into study folder.
    mech_template_setup(study_id)

    # 4. Run parameter estimation:
    #C1_init = 2.0   # Default value.

    if a_p_toggle == 1:
        # Passive parameter estimation only.
        #[c1_opt, opt_val] = optimise_passive_main(study_id, study_frame, pressure, data_idx, c1_init,
        #                                          forward_solve_toggle)
        # Evaluate displacement by projecting surface data onto high stiffness model.
        #optimised_passive_displacement_mse(study_id, study_frame)
        # Perform passive parameter sweep
        #time.sleep(5)
        #setup_log_file(study_id, 'Output_Passive_Parameter_Sweep.log')
        #optimise_passive_parameter_sweep_sensible_projection(study_id, study_frame)
        optimise_passive_parameter_sweep_ed_only(study_id, study_frame)
        #return c1_opt, opt_val, pressure[0]
    elif a_p_toggle == 2:
        # apex_acti = optimise_determine_apical_activation(study_id, study_frame, pressure, data_idx)
        apex_acti = 1
        print 'LOG_STATUS: Apical Activation switch is: ', str(apex_acti)
        # Active parameter estimation only.
        optimise_active_main(study_id, study_frame, pressure, data_idx, apex_acti)
        # 5. Evaluate TCa sensitivity at each optimised frame.
        active_sensitivity_evaluate(study_id, study_frame)
    elif a_p_toggle == 3:
        # Passive and active parameter estimations.
        optimise_passive_main(study_id, study_frame, pressure, data_idx, c1_init, forward_solve_toggle)
        apex_acti = optimise_determine_apical_activation(study_id, study_frame, pressure, data_idx)
        optimise_active_main(study_id, study_frame, pressure, data_idx, apex_acti)
        active_sensitivity_evaluate(study_id, study_frame)


def main():
    # Main commands begin here.
    os.system('rm *.*~')  # Remove all temporary files in Main folder.

    # Top level control of the entire simulation:
    # Command line inputs:
    idx = int(sys.argv[1])  # Study number
    debug_toggle = int(sys.argv[2])  # Debug toggle: 1: write output to text file, 0: write output to terminal.
    forward_toggle = int(sys.argv[3])  # Forward solve toggle: 1 - do passive initial solve, 0 - don't do it.
    a_p_toggle = int(sys.argv[4])  # Run: 1. Passive only, 2. Active only, 3. Passive & active.

    # Extract the important frame numbers for all studies
    [study_id, study_frame, c1_init] = setup_get_frames(idx)

    # Output to screen
    print ''
    print '********************************************************'
    print '          Analysing study ' + study_id
    print '********************************************************'
    print ''

    # Run main program.
    #main_program(study_id, study_frame, c1_init, debug_toggle, forward_toggle, a_p_toggle)
    [c1_opt, opt_val, edp] = main_program(study_id, study_frame, c1_init, debug_toggle, forward_toggle, a_p_toggle)
    #[eps, d2Error, num_frames] = hessian_test_epsilon(c1_opt, opt_val, study_id, study_frame)
    #results_write_to_text_file(idx, c1_opt, opt_val, eps, d2Error, num_frames, edp)

if __name__ == '__main__':
    main()
