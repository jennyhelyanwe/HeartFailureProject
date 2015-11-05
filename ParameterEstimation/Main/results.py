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
    f_w.write('Pressure (kPa)\tVolume (mm^3)\tMSE Endo (mm^2)\tMSE Epi (mm^2)\tMSE total (mm^2)\n')
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

    f_w = open(os.environ['RESULTS']+study_id+'/Passive_Regional_Errors_'+name+'.txt', 'w+')
    f_w.write('Epi_apex\tEpi_midlow\tEpi_mid\tEpi_base\tEndo_apex\tEndo_midlow\tEndo_mid\tEndo_base\n')
    f_w.write('0\t0\t0\t0\t0\t0\t0\t0\n')
    for i in range(1, len(idx)):
        fr = open(status+'Error/EpiError_'+str(idx[i])+'.exdata', 'r')
        data = fr.readlines()
        j = 9
        epi_apex = []
        while j < 308:
            j += 3
            temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
            epi_apex.append(numpy.linalg.norm(temp))
        epi_apex_ave = numpy.mean(epi_apex)
        epi_midlow = []
        while j < 1508:
            j += 3
            temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
            epi_midlow.append(numpy.linalg.norm(temp))
        epi_midlow_ave = numpy.mean(epi_midlow)
        epi_mid = []
        while j < 2708:
            j += 3
            temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
            epi_mid.append(numpy.linalg.norm(temp))
        epi_mid_ave = numpy.mean(epi_mid)
        epi_base = []
        while j < 3908:
            j += 3
            temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
            epi_base.append(numpy.linalg.norm(temp))
        epi_base_ave = numpy.mean(epi_base)

        fr = open(status+'Error/EndoError_'+str(idx[i])+'.exdata', 'r')
        data = fr.readlines()
        j = 9
        endo_apex = []
        while j < 308:
            j += 3
            temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
            endo_apex.append(numpy.linalg.norm(temp))
        endo_apex_ave = numpy.mean(endo_apex)
        endo_midlow = []
        while j < 1508:
            j += 3
            temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
            endo_midlow.append(numpy.linalg.norm(temp))
        endo_midlow_ave = numpy.mean(endo_midlow)
        endo_mid = []
        while j < 2708:
            j += 3
            temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
            endo_mid.append(numpy.linalg.norm(temp))
        endo_mid_ave = numpy.mean(endo_mid)
        endo_base = []
        while j < 3908:
            j += 3
            temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
            endo_base.append(numpy.linalg.norm(temp))
        endo_base_ave = numpy.mean(endo_base)
        f_w.write(str(epi_apex_ave)+'\t'+str(epi_midlow_ave)+'\t'+str(epi_mid_ave)+'\t'+str(epi_base_ave)+'\t'+
                  str(endo_apex_ave)+'\t'+str(endo_midlow_ave)+'\t'+str(endo_mid_ave)+'\t'+str(endo_base_ave)+'\n')

    if name == 'Opt':
        copy('LV_CubicGuc.ipmate', os.environ['RESULTS']+study_id+'/LV_CubicGuc_Opt.ipmate')

    """
    for i in range(1, len(idx)):
        with open(os.environ['RESULTS']+study_id+'/Stress_ER_'+str(idx[i])+'.txt', 'w') as f_w:
            f_w.write('Fibre extension ratio: \n')
            filename = 'OptimisedStressStrain/ER_'+str(idx[i])+'.exdata'
            [ER_average] = average_stress_strain_elem(filename, component)
            f_w.write('%f', ER_average)
    """



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
            TCa = data[17].split()[-1]

            f_w.write(str(p)+'\t'+str(TCa)+'\t'+str(v)+'\t'+str(e_endo)+'\t'+str(e_epi)+'\t'+str(e_tot)+'\n')

    with open(os.environ['RESULTS']+study_id+'/Active_Regional_Errors_'+name+'.txt', 'w+') as f_w:
        f_w.write('Epi_apex\tEpi_midlow\tEpi_mid\tEpi_base\tEndo_apex\tEndo_midlow\tEndo_mid\tEndo_base\n')
        for i in range(1, len(idx)):
            fr = open(status+'Error/EpiError_'+str(idx[i])+'.exdata', 'r')
            data = fr.readlines()
            j = 9
            epi_apex = []
            while j < 308:
                j += 3
                temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
                epi_apex.append(numpy.linalg.norm(temp))
            epi_apex_ave = numpy.mean(epi_apex)
            epi_midlow = []
            while j < 1508:
                j += 3
                temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
                epi_midlow.append(numpy.linalg.norm(temp))
            epi_midlow_ave = numpy.mean(epi_midlow)
            epi_mid = []
            while j < 2708:
                j += 3
                temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
                epi_mid.append(numpy.linalg.norm(temp))
            epi_mid_ave = numpy.mean(epi_mid)
            epi_base = []
            while j < 3908:
                j += 3
                temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
                epi_base.append(numpy.linalg.norm(temp))
            epi_base_ave = numpy.mean(epi_base)

            fr = open(status+'Error/EndoError_'+str(idx[i])+'.exdata', 'r')
            data = fr.readlines()
            j = 9
            endo_apex = []
            while j < 308:
                j += 3
                temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
                endo_apex.append(numpy.linalg.norm(temp))
            endo_apex_ave = numpy.mean(endo_apex)
            endo_midlow = []
            while j < 1508:
                j += 3
                temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
                endo_midlow.append(numpy.linalg.norm(temp))
            endo_midlow_ave = numpy.mean(endo_midlow)
            endo_mid = []
            while j < 2708:
                j += 3
                temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
                endo_mid.append(numpy.linalg.norm(temp))
            endo_mid_ave = numpy.mean(endo_mid)
            endo_base = []
            while j < 3908:
                j += 3
                temp = [float(data[j].split()[0]), float(data[j].split()[1]), float(data[j].split()[2])]
                endo_base.append(numpy.linalg.norm(temp))
            endo_base_ave = numpy.mean(endo_base)
            f_w.write(str(epi_apex_ave)+'\t'+str(epi_midlow_ave)+'\t'+str(epi_mid_ave)+'\t'+str(epi_base_ave)+'\t'+
                      str(endo_apex_ave)+'\t'+str(endo_midlow_ave)+'\t'+str(endo_mid_ave)+'\t'+str(endo_base_ave)+'\n')
#
#=======================================================================================================================
#

