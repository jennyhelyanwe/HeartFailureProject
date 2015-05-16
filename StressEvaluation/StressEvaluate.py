# This script evaluates the equatorial maximum principal stress at the maximum activation point in the cardiac cycle
# It then generates a 2D figure with the first layer of gauss points from the equatorial element. 

import os, sys
import numpy 
from numpy import array
from readFiles import *

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

# Get TCa transient
if os.path.exists('ActivationOpt/'):
    os.chdir('ActivationOpt/')
    a_folder = 'ActivationOpt/'
elif os.path.exists('OptimisedActivation/'):
    os.chdir('OptimisedActivation/')
    a_folder = 'OptimisedActivation/'
all_files = os.listdir('./')
TCa = numpy.zeros(len(all_files))
for f in all_files:
    fname = f.split('.')[0]
    frame = int(fname.split('_')[1])
    fid = open(f, 'r')
    temp = fid.readlines()
    if temp[0].split()[2] == '2.1':
        toggle = True
        TCa[frame-1] = float(temp[15].split()[6])
    else:
        toggle = False
        TCa[frame-1] = float(temp[21].split()[4])
        # If the old version, then we need to update correct the ipacti file.
        filename= 'TCa_'+str(frame)+'.ipacti'
        write_ipacti_file(TCa[frame-1], filename)

# Get frame number for maximum TCa
idx = numpy.argmax(TCa)
max_TCa = numpy.max(TCa)
print max_TCa

# Edit command file to read the correct frame in.
with open('../StressStrainEvaluate.com', 'r') as file:
    data = file.readlines()

#if toggle:
#    data[15] = 'fem def acti;r;ActivationOpt/TCa_'+str(idx+1)+'\n'
#else:
#    data[15] = 'fem def acti;r;ActivationOpt/calcium0xx_'+str(idx+1)+'\n'
data[15] = 'fem def acti;r;'+a_folder+'TCa_'+str(idx+1)+'\n'
if a_folder == 'ActivationOpt/':
    data[6] = 'fem def node;r;CAP_DS_humanLV\n'
    data[7] = 'fem def elem;r;CAP_DS_humanLV\n'
    data[16] = 'fem def init;r;PostContraction_BC/CurrentContracted_'+str(idx+1)+'\n'
else:
    data[6] = 'fem def node;r;DSWall\n'
    data[7] = 'fem def elem;r;DSWall\n'
    data[16] = 'fem def init;r;WarmStartSolution/CurrentContracted_'+str(idx+1)+'\n'

data[20] = 'fem list stress;FullStress_'+str(idx+1)+' full\n'
data[21] = 'fem list stress;PassiveStress_'+str(idx+1)+' passive full\n'
data[22] = 'fem list stress;ActiveStress_'+str(idx+1)+' active full\n'
data[23] = 'fem list strain;Strain_'+str(idx+1)+' full\n'

with open('../StressStrainEvaluate.com', 'w') as file:
    file.writelines(data)

os.chdir('../')
os.system('cm StressStrainEvaluate.com')
#################################################################################
# Read in total principal stress at all gauss points.
filename = 'FullStress_'+str(idx+1)+'.opstre'
[x, xi, t_mps, t_fs_2PK, t_fs_C] = readOpstre(filename)

## Write out mps
equatorial_xi = xi[8:12]
equatorial_x = x[8:12]
equatorial_t_mps = t_mps[8:12]
TOL = 1e-3
with open('TotalStress2D_'+str(idx+1)+'.exdata', 'w') as file:
    file.write(' Group name: Total_mps\n')
    file.write(' #Fields=2\n')
    file.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
    file.write('   x.  Value index= 1, #Derivatives=0\n')
    file.write('   y.  Value index= 2, #Derivatives=0\n')
    file.write('   z.  Value index= 3, #Derivatives=0\n')
    file.write(' 1) total_mps, coordinate, rectangular cartesian, #Components=1\n')
    file.write('   1.  Value index= 1, #Derivatives=0\n')
    k = 0
    eq_t_col = []
    for i in range(0, len(equatorial_xi)):
        for j in range(0,64):
            if abs(equatorial_xi[i][j][1]-0.93057) < TOL:
                # Write this point to exdata file.
                k = k + 1
                file.write(' Node:\t'+str(k)+'\n')
                file.write(' '+str(equatorial_x[i][j][0])+'\t'+str(equatorial_x[i][j][1])+'\t'+str(equatorial_x[i][j][2])+'\n')
                file.write(' '+str(equatorial_t_mps[i][j])+'\n')
            eq_t_col.append(float(equatorial_t_mps[i][j]))

avg_mps_eq = numpy.mean(eq_t_col)
print 'Total mps in equatorial elements:'
print avg_mps_eq
std_mps_eq = numpy.std(eq_t_col)
print 'std:'
print std_mps_eq

t_col = []
for i in range (0, len(t_mps)):
    for j in range(0, 64):
        t_col.append(float(t_mps[i][j]))
avg_mps = numpy.mean(t_col)
std_mps = numpy.std(t_col)
print 'Total mps over entire model'
print avg_mps
print 'std: '
print std_mps

## Write out fibre stress
eq_fs_C_col = []
t_fs_C_col = []
for i in range(0, len(t_fs_C)):
    for j in range(0, 64):
        t_fs_C_col.append(float(t_fs_C[i][j]))
        if (i >= 8) & (i < 12):
            eq_fs_C_col.append(float(t_fs_C[i][j]))

avg_fs_C_eq = numpy.mean(eq_fs_C_col)
std_fs_C_eq = numpy.std(eq_fs_C_col)
print 'Total fibre stress in equatorial elements'
print avg_fs_C_eq
print 'std: '
print std_fs_C_eq

avg_fs_C = numpy.mean(t_fs_C_col)
std_fs_C = numpy.std(t_fs_C_col)
print 'Total fibre stress over entire model'
print avg_fs_C
print 'std: '
print std_fs_C

# Run CMGUI command file to save picture of stress.
with open('EquatorialVisualise.com', 'r') as file:
    data = file.readlines()

if a_folder == 'ActivationOpt/':
    data[4] = 'gfx read node output_VisualStep_Opt/Contracted_'+str(idx+1)+' region model\n'
    data[5] = 'gfx read elem output_VisualStep_Opt/Contracted_'+str(idx+1)+' region model\n'
else:
    data[3] = 'gfx read node OptimisedExfile/LVInflation_'+str(idx+1)+' region model\n'
    data[4] = 'gfx read elem OptimisedExfile/LVInflation_'+str(idx+1)+' region model\n'
data[5] = 'gfx read data TotalStress2D_'+str(idx+1)+'\n'
data[29] = 'gfx print png file TotalStress2D.png window 1;\n'

with open('EquatorialVisualise.com', 'w') as file:
    file.writelines(data)

os.system('cmgui EquatorialVisualise.com')

#################################################################################
# Read in total principal stress at all gauss points.
filename = 'ActiveStress_'+str(idx+1)+'.opstre'
[x, xi, a_mps, a_fs_2PK, a_fs_C] = readActiveOpstre(filename)

## Write out maximum principal stresses.
equatorial_xi = xi[8:12]
equatorial_x = x[8:12]
equatorial_a_mps = a_mps[8:12]
TOL = 1e-3
with open('ActiveStress2D_'+str(idx+1)+'.exdata', 'w') as file:
    file.write(' Group name: Active_mps\n')
    file.write(' #Fields=2\n')
    file.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
    file.write('   x.  Value index= 1, #Derivatives=0\n')
    file.write('   y.  Value index= 2, #Derivatives=0\n')
    file.write('   z.  Value index= 3, #Derivatives=0\n')
    file.write(' 1) total_mps, coordinate, rectangular cartesian, #Components=1\n')
    file.write('   1.  Value index= 1, #Derivatives=0\n')
    k = 0
    eq_a_col = []
    for i in range(0, len(equatorial_xi)):
        for j in range(0,64):
            if abs(equatorial_xi[i][j][1]-0.93057) < TOL:
                # Write this point to exdata file.
                k = k + 1
                file.write(' Node:\t'+str(k)+'\n')
                file.write(' '+str(equatorial_x[i][j][0])+'\t'+str(equatorial_x[i][j][1])+'\t'+str(equatorial_x[i][j][2])+'\n')
                file.write(' '+str(equatorial_a_mps[i][j])+'\n')
            eq_a_col.append(float(equatorial_a_mps[i][j]))

avg_a_mps_eq = numpy.mean(eq_a_col)
std_a_mps_eq = numpy.std(eq_a_col)
print 'Active mps in equatorial elements'
print avg_a_mps_eq
print 'std'
print std_a_mps_eq

a_col = []
for i in range(0, len(a_mps)):
    for j in range(0, 64):
        a_col.append(float(a_mps[i][j]))

avg_a_mps = numpy.mean(a_col)
std_a_mps = numpy.std(a_col)
print 'Active mps over entire model'
print avg_a_mps
print 'std: '
print std_a_mps

## Write out fibre stresses.
a_fs_C_col = []
eq_a_fs_C_col = []
for i in range(0, len(a_fs_C)):
    for j in range(0, 64):
        a_fs_C_col.append(float(a_mps[i][j]))
        if (i>=8) & (i <12):
            eq_a_fs_C_col.append(float(a_fs_C[i][j]))

eq_avg_a_fs_C = numpy.mean(eq_a_fs_C_col)
eq_std_a_fs_C = numpy.std(eq_a_fs_C_col)
print 'Active fibre stress in equatorial elements:'
print eq_avg_a_fs_C
print 'std:'
print eq_std_a_fs_C


avg_a_fs_C = numpy.mean(a_fs_C_col)
std_a_fs_C = numpy.std(a_fs_C_col)
print 'Active fibre stress over entire model:'
print avg_a_fs_C
print 'std:'
print std_a_fs_C



# Run CMGUI command file to save picture of stress.
with open('EquatorialVisualise_Active.com', 'r') as file:
    data = file.readlines()
if a_folder == 'ActivationOpt/':
    data[4] = 'gfx read node output_VisualStep_Opt/Contracted_'+str(idx+1)+' region model\n'
    data[5] = 'gfx read elem output_VisualStep_Opt/Contracted_'+str(idx+1)+' region model\n'
else:
    data[3] = 'gfx read node OptimisedExfile/LVInflation_'+str(idx+1)+' region model\n'
    data[4] = 'gfx read elem OptimisedExfile/LVInflation_'+str(idx+1)+' region model\n'
data[5] = 'gfx read data ActiveStress2D_'+str(idx+1)+'\n'
data[29] = 'gfx print png file ActiveStress2D.png window 1;\n'

with open('EquatorialVisualise_Active.com', 'w') as file:
    file.writelines(data)

os.system('cmgui EquatorialVisualise_Active.com')
#################################################################################
# Read in total principal stress at all gauss points.
filename = 'PassiveStress_'+str(idx+1)+'.opstre'
[x, xi, p_mps, p_fs_2PK, p_fs_C] = readOpstre(filename)

# Pick out the gauss points in elements 9 to 12 which have a xi2 value of 0.93057
equatorial_xi = xi[8:12]
equatorial_x = x[8:12]
equatorial_p_mps = p_mps[8:12]
TOL = 1e-3
with open('PassiveStress2D_'+str(idx+1)+'.exdata', 'w') as file:
    file.write(' Group name: Passive_mps\n')
    file.write(' #Fields=2\n')
    file.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
    file.write('   x.  Value index= 1, #Derivatives=0\n')
    file.write('   y.  Value index= 2, #Derivatives=0\n')
    file.write('   z.  Value index= 3, #Derivatives=0\n')
    file.write(' 1) total_mps, coordinate, rectangular cartesian, #Components=1\n')
    file.write('   1.  Value index= 1, #Derivatives=0\n')
    k = 0
    eq_p_col = []
    for i in range(0, len(equatorial_xi)):
        for j in range(0,64):
            if abs(equatorial_xi[i][j][1]-0.93057) < TOL:
                # Write this point to exdata file.
                k = k + 1
                file.write(' Node:\t'+str(k)+'\n')
                file.write(' '+str(equatorial_x[i][j][0])+'\t'+str(equatorial_x[i][j][1])+'\t'+str(equatorial_x[i][j][2])+'\n')
                file.write(' '+str(equatorial_p_mps[i][j])+'\n')
                eq_p_col.append(float(equatorial_p_mps[i][j]))

avg_p_mps_eq = numpy.mean(eq_p_col)
std_p_mps_eq = numpy.std(eq_p_col)
print 'Passive mps in equatorial elements'
print avg_p_mps_eq
print 'std'
print std_p_mps_eq

p_col = []
for i in range(0, len(p_mps)):
    for j in range(0, 64):
        p_col.append(float(p_mps[i][j]))

avg_p_mps = numpy.mean(p_col)
std_p_mps = numpy.std(p_col)
print 'Passive mps over entire model'
print avg_p_mps
print 'std: '
print std_p_mps

## Write out fibre stresses.
p_fs_C_col = []
eq_p_fs_C_col = []
for i in range (0, len(p_fs_C)):
    for j in range(0, 64):
        p_fs_C_col.append(float(p_fs_C[i][j]))
        if (i >=8) & (i < 12):
            eq_p_fs_C_col.append(float(p_fs_C[i][j]))

eq_avg_p_fs_C = numpy.mean(eq_p_fs_C_col)
eq_std_p_fs_C = numpy.std(eq_p_fs_C_col)
print 'Passive fibre stress in equatorial elements:'
print eq_avg_p_fs_C
print 'std:'
print eq_std_p_fs_C


avg_p_fs_C = numpy.mean(p_fs_C_col)
std_p_fs_C = numpy.std(p_fs_C_col)
print 'Passive fibre stress over entire model:'
print avg_p_fs_C
print 'std:'
print std_p_fs_C

# Run CMGUI command file to save picture of stress.
with open('EquatorialVisualise_Passive.com', 'r') as file:
    data = file.readlines()

if a_folder == 'ActivationOpt/':
    data[4] = 'gfx read node output_VisualStep_Opt/Contracted_'+str(idx+1)+' region model\n'
    data[5] = 'gfx read elem output_VisualStep_Opt/Contracted_'+str(idx+1)+' region model\n'
else:
    data[3] = 'gfx read node OptimisedExfile/LVInflation_'+str(idx+1)+' region model\n'
    data[4] = 'gfx read elem OptimisedExfile/LVInflation_'+str(idx+1)+' region model\n'
data[6] = 'gfx read data PassiveStress2D_'+str(idx+1)+'\n'
data[29] = 'gfx print png file PassiveStress2D.png window 1;\n'

with open('EquatorialVisualise_Passive.com', 'w') as file:
    file.writelines(data)

os.system('cmgui EquatorialVisualise_Passive.com')

#################################################################################
# Read in strains at all gauss points.
filename = 'Strain_'+str(idx+1)+'.opstra'
[x, xi, exr] = readOpstra(filename)

exr_col = []
eq_exr_col = []
for i in range(0, len(exr)):
    for j in range(0, 64):
        exr_col.append(float(exr[i][j]))
        if (i >= 8) & (i < 12):
            eq_exr_col.append(float(exr[i][j]))

avg_exr = numpy.mean(exr_col)
std_exr = numpy.std(exr_col)
print 'Fibre extension ratio over entire model'
print avg_exr
print 'std:'
print std_exr

eq_avg_exr = numpy.mean(eq_exr_col)
eq_std_exr = numpy.std(eq_exr_col)
print 'Fibre extension ratio in equatorial elements'
print eq_avg_exr
print 'std:'
print eq_std_exr


