__author__ = 'zwan145'

import os
from passive_mechanics import *
from active_mechanics import *


def results_passive_generate(study_id, study_frame, name, status):
    # Copy optimised parameters to results folder
    os.chdir(os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/PassiveMechanics')

    # Generate PV curve for optimised models with fitting error at each frame.
    f_w = open(os.environ['RESULTS']+study_id+'/Passive_PV_Errors_'+name+'.txt', 'w+')

    # Write DS volume and error of zero with DS pressure
    f_r = open('output_cavity_volume/LVCavityInit.opelem', 'r')
    data = f_r. readlines()
    v = data[5].split()[3]
    f_r.close()
    #f_w.write('Pressure (kPa)\tVolume (mm^3)\tMSE Endo (mm^2)\tMSE Epi (mm^2)\tMSE total (mm^2)\n')
    f_w.write('0.0\t'+str(v)+'\t0.0\t0.0\t0.0\n')

    idx = passive_loop_index(study_frame)
    for i in range(1, len(idx)):
        f_r = open('pressure/pressure_'+str(idx[i])+'.opvari','r')
        data = f_r.readlines()
        p = data[1006].split()[3]
        f_r.close()

        f_r = open(status+'CavityVolume/LVCavity_'+str(idx[i])+'.opelem', 'r')
        data = f_r.readlines()
        v = data[5].split()[3]
        f_r.close()

        f_r = open(status+'Error/EndoError_'+str(idx[i])+'.opdata', 'r')
        data = f_r.readlines()
        e_endo = data[7].split()[5]
        e_endo = float(e_endo)**2
        f_r.close()

        f_r = open(status+'Error/EpiError_'+str(idx[i])+'.opdata', 'r')
        data = f_r.readlines()
        e_epi = data[7].split()[5]
        e_epi = float(e_epi)**2
        f_r.close()

        e_tot = (e_endo + e_epi)/2

        f_w.write(str(p)+'\t'+str(v)+'\t'+str(e_endo)+'\t'+str(e_epi)+'\t'+str(e_tot)+'\n')

    f_w.close()

    if name == 'Opt':
        copy('LV_CubicGuc.ipmate', os.environ['RESULTS']+study_id+'/LV_CubicGuc_Opt.ipmate')
#
#=======================================================================================================================
#


def results_active_generate(study_id, study_frame, name, status):
    # Copy optimised parameters to results folder
    os.chdir(os.environ['STUDIES']+study_id+'/LVMechanics'+study_id+'/ActiveMechanics')

    with open(os.environ['RESULTS']+study_id+'/Active_PV_Errors_'+name+'.txt', 'w+') as f_w:
        f_w.write('Pressure (kPa)\tTCa (kPa)\tVolume (mm^3)\tEndo error (mm^2)\tEpi error (mm^2)\tTotal error (mm^2)\n')
        idx = active_loop_index(study_frame)
        for i in range(1, len(idx)):
            f_r = open('pressure/pressure_'+str(idx[i])+'.opvari','r')
            data = f_r.readlines()
            p = data[1006].split()[3]
            f_r.close()

            f_r = open(status+'CavityVolume/LVCavity_'+str(idx[i])+'.opelem', 'r')
            data = f_r.readlines()
            v = data[5].split()[3]
            f_r.close()

            f_r = open(status+'Error/EndoError_'+str(idx[i])+'.opdata', 'r')
            data = f_r.readlines()
            e_endo = data[7].split()[5]
            e_endo = float(e_endo)**2
            f_r.close()

            f_r = open(status+'Error/EpiError_'+str(idx[i])+'.opdata', 'r')
            data = f_r.readlines()
            e_epi = data[7].split()[5]
            e_epi = float(e_epi)**2
            f_r.close()

            e_tot = (e_endo + e_epi)/2

            f_r = open('OptimisedActivation/TCa_'+str(idx[i])+'.ipacti', 'r')
            data = f_r.readlines()
            TCa = data[15].split()[6]

            f_w.write(str(p)+'\t'+str(TCa)+'\t'+str(v)+'\t'+str(e_endo)+'\t'+str(e_epi)+'\t'+str(e_tot)+'\n')
#
#=======================================================================================================================
#

