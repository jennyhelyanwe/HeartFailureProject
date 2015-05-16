import os, sys
import numpy
from readFiles import *

"""
Evaluate T_a transient
----------------------
Author: Z.J. Wang
Date: 07 April 2015
----------------------

This script controls the CMISS command file which will read in each frame
in the cardiac cycle and evaluate the active fibre stress averaged over the
entire model as well as the standard deviation. It then writes this into a
file ready to be imported into Microsoft Excel to be plotted.

"""

def write_ipacti_file(TCa, filename):
    with open(filename, 'w') as file:
        file.write(' CMISS Version 2.1  ipacti File Version 2\n')
        file.write(' Heading:\n')
        file.write('\n')
        file.write(' Specify type of contraction model [2]:\n')
        file.write('   (1) SS tension-length-Ca relation (set Cai)\n')
        file.write('   (2) Get active stress from cell/grid model  (e.g. Hunter/McCulloch/ter Keurs model)\n')
        file.write('   (3) Enhance 2nd PKST by user-defined functions\n')
        file.write('   (4) SS TCa-length relation (set TCa)\n')
        file.write('    4\n')
        file.write(' Enter non-dimensional slope parameter (beta) [1.45]: 0.14500E+01\n')
        file.write(' Specify whether the initial tension level TCa is [1]:\n')
        file.write('   (1) Constant spatially\n')
        file.write('   (2) Piecewise constant (defined by elements)\n')
        file.write('   (3) Defined by Gauss points\n')
        file.write('    1\n')
        file.write(' Enter initial tension level TCa [0]: '+str(TCa)+'\n') 

# Get T_Ca transient
if os.path.exists('ActivationOpt/'):
    os.chdir('ActivationOpt/')
    a_folder = 'ActivationOpt/'
    os.system('rm TCa_*')
elif os.path.exists('OptimisedActivation/'):
    os.chdir('OptimisedActivation/')
    a_folder = 'OptimisedActivation/'
all_files = os.listdir('./')
TCa = []
for f in all_files:
    fname = f.split('.')[0]
    frame = int(fname.split('_')[1])
    fid = open(f, 'r')
    temp = fid.readlines()
    if temp[0].split()[2] == '2.1':
        toggle = True
        TCa.append(float(temp[15].split()[6]))
    else:
        toggle = False
        if fname.split('_')[0] == 'calcium0xx':
            TCa.append(float(temp[21].split()[4]))
            # If the old version, then we need to update correct the ipacti file.
            filename= 'TCa_'+str(frame)+'.ipacti'
            write_ipacti_file(float(temp[21].split()[4]), filename)

# For each frame, read in solution and TCa, evaluate stress.
os.chdir('../')
with open('Ta.txt', 'w') as f_w:
    # Write labels for values.
    f_w.write('Frame number\tTot. fib.\tstd\tPass. fib.\tstd\tAct. fib. \tstd\tFib. ER\tstd\n')
    for i in range(1, len(TCa)):  # Start from frame 2.
        # Total fibre stress
        [x, xi, t_mps, fs_2PK, fs_C] = readOpstre('OptimisedStressStrain/total_stress_'+str(i+1)+'.opstre')
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
