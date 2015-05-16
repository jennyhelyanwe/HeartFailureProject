import os, sys
import numpy
from readFiles import *

"""
Evaluate stress transient
----------------------
Author: Z.J. Wang
Date: 07 April 2015
----------------------

This script controls the CMISS command file which will read in each frame
in the cardiac cycle and evaluate the total, passive and active fibre stress
averaged over the entire model as well as the standard deviation. It then
writes this into a file ready to be imported into Microsoft Excel.

"""

"""
Passive stress and strain transients
===================================
"""

os.chdir('WarmStartSolution')
all_files = os.listdir('./')
idx = []
for f in all_files:
    fname = f.split('.')[0]
    frame = int(fname.split('_')[1])
    idx.append(frame)
idx = sorted(idx)
# For each frame, read in solution and TCa, evaluate stress.
os.chdir('../')
with open('StressStrainTransients.txt', 'w') as f_w:
    # Write labels for values. 
    f_w.write('Frame number\tTot. fib.\tstd\tPass. fib.\tstd\tAct. fib. \tstd\tFib. ER\tstd\n')
    for i in range(0, len(idx)):  # Start from frame 2.
        # Total fibre stress
        [x, xi, t_mps, fs_2PK, fs_C] = readOpstre('OptimisedStressStrain/total_stress_'+str(idx[i])+'.opstre')
        fs = []
        for j in range(0, len(fs_C)):
            for k in range(0, 64):
                fs.append(float(fs_C[j][k]))
        avg = numpy.mean(fs)
        std = numpy.std(fs)

        # Active fibre stress
        [x, xi, a_t_mps, a_fs_2PK, a_fs_C] = readOpstre('OptimisedStressStrain/active_stress_'+str(i+1)+'.opstre')
        fs_a = []
        for j in range(0, len(a_fs_C)):
            for k in range(0, 64):
                fs_a.append(float(a_fs_C[j][k]))
        avg_a = numpy.mean(fs_a)
        std_a = numpy.std(fs_a)

        # Passive fibre stress
        [x, xi, p_t_mps, p_fs_2PK, p_fs_C] = readOpstre('OptimisedStressStrain/passive_stress_'+str(i+1)+'.opstre')
        fs_p = []
        for j in range(0, len(p_fs_C)):
            for k in range(0, 64):
                fs_p.append(float(p_fs_C[j][k]))
        avg_p = numpy.mean(fs_p)
        std_p = numpy.std(fs_p)

        # Read in strains at all gauss points.
        filename = 'Strain_'+str(i+1)+'.opstra'
        [x, xi, exr] = readOpstra(filename)
        vector = []
        for i in range(0, len(exr)):
            for j in range(0, 64):
                vector.append(float(exr[i][j]))
        avg_exr = numpy.mean(vector)
        std_exr = numpy.std(vector)

        # Write to file.
        frame_no = i+1
        f_w.write('%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.3f\t%.3f\n' % (frame_no, avg, std, avg_a, std_a, avg_p, std_p, avg_exr, std_exr))


        """
        with open('ActiveStress.com', 'r') as file:
            data=file.readlines()
        data[15] = 'fem def acti;r;'+a_folder+'TCa_'+str(i+1)+'\n'
        print a_folder
        if a_folder == 'ActivationOpt/':
            data[6] = 'fem def node;r;CAP_DS_humanLV\n'
            data[7] = 'fem def elem;r;CAP_DS_humanLV\n'
            data[16] = 'fem def init;r;PostContraction_BC/CurrentContracted_'+str(i+1)+'\n'
        else:
            data[6] = 'fem def node;r;DSWall\n'
            data[7] = 'fem def elem;r;DSWall\n'
            data[16] = 'fem def init;r;WarmStartSolution/CurrentContracted_'+str(i+1)+'\n'
        data[20] = 'fem list stress;ActiveStress_'+str(i+1)+' active full\n'

        with open('ActiveStress.com', 'w') as file:
            file.writelines(data)

        os.system('cm ActiveStress.com')

        # Read active stress of current file from opstre file.
        filename = 'ActiveStress_'+str(i+1)+'.opstre'
        [x, xi, a_mps, a_fs_2PK, a_fs_C] = readActiveOpstre(filename)
        fs_vector = []
        for j in range(0, len(a_fs_C)):
            for k in range(0, 64):
                fs_vector.append(float(a_fs_C[j][k]))

        avg = numpy.mean(fs_vector)
        std = numpy.std(fs_vector)
        frame_no = i+1
        """
