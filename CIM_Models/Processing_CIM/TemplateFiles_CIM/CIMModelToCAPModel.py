import sys
import os
from numpy.ma import zeros
from numpy.numarray.numerictypes import Float
from numpy.numarray.numerictypes import Int

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="jchu014"
__date__ ="$Mar 18, 2010 12:15:45 PM$"

def read_map(filename):
    f = open(filename)

    i = 0
    while i < 2:
        line = f.readline()
        print line
        i = i + 1
    
    line = f.readline() # RUA
    dimensions = line.split()
    numRow = int(dimensions[1])
    numCol = int(dimensions[2])
    numNonZero = int(dimensions[3])
    
    print 'numRow = ', numRow,'numCol = ',  numCol\
    , 'numNonZero = ', numNonZero

    map = zeros([numCol,numRow], Float)

    line = f.readline()
    print line

    columnIndex = []
    while len(columnIndex) <= numCol:
        line = f.readline()
        temp = line.split()
        columnIndex.extend(temp)
#    print columnIndex

    rowIndex = []
    while len(rowIndex) < numNonZero:
        line = f.readline()
        temp = line.split()
        rowIndex.extend(temp)
#    print rowIndex

    valuesString = f.read()
    values = valuesString.split()
    print "num values", len(values)

    count = 0
    for i in range(numCol):
        endIndex = int(columnIndex[i+1])-1
        indexIntoRowIndex = int(columnIndex[i])-1
#        print "endIndex", endIndex
        while indexIntoRowIndex < endIndex:
#            print indexIntoRowIndex,
#            print count,
            row = int(rowIndex[indexIntoRowIndex])-1
#            print row
            map[i][row] = float(values[count])
            count = count + 1
            indexIntoRowIndex = indexIntoRowIndex + 1
    return map

def convert_to_hermite(filename, map0, map1, map2):
    f = open(filename) # CIM model file
    line = f.readline() #'FocalLength'
    print line
    focus = float(f.readline())

    for i in range(10):
        line = f.readline()
        print i, line

    globalParam_0 = [] # lambda
    for i in range(134):
        line = f.readline()
        print i, line
        globalParam_0.append(float(line))

    for i in range(4):
        line = f.readline()
        print i, line

    globalParam_1 = [] # mu
    for i in range(40):
        globalParam_1.append(float(f.readline()))


    for i in range(4):
        line = f.readline()
        print i, line

    globalParam_2 = [] # theta
    for i in range(40):
        globalParam_2.append(float(f.readline()))
    print globalParam_2

    # Multiply

    localParam_0 = zeros([512], Float)
    for r in range(512):
        temp = 0.0
        for c in range(134):
            temp += map0[c][r] * globalParam_0[c]
        localParam_0[r] = temp

    localParam_1 = zeros([128], Float)
    for r in range(128):
        temp = 0.0
        for c in range(40):
            temp += map1[c][r] * globalParam_1[c]
        localParam_1[r] = temp

    localParam_2 = zeros([128], Float)
    for r in range(128):
        temp = 0.0
        for c in range(40):
            temp += map2[c][r] * globalParam_2[c]
        localParam_2[r] = temp

    print localParam_2

    #Now print the global node parameters
    #
    #first rearrange the vector of local parameters to get the global parameters (for Hermite basis)
    #this could altenatively done by multipling the local parameter vector with the inverse of Global_to_local map
    # (which can be found from CimModelLVPS4x4_1.map)
    # but here we just use the map to extract the info to map local param back to global map as follows

    index = zeros([40], Int)
    for c in range(40):
        for r in range(160):
            if map2[c][r] == 1.0:
                index[c] = r
                print c + 1, ": ", r + 1
                break
    #now print

    outfile = filename + ".exnode"
    of = open(outfile, 'w')
    of.write("Region: /heart\n")
    of.write("#Fields=1\n")
    line = "1) coordinates, coordinate, prolate spheroidal, focus= %s, #Components=3\n" % focus
    of.write(line)
    of.write("  lambda.  Value index= 1, #Derivatives= 3 (d/ds1,d/ds2,d2/ds1ds2)\n")
    of.write("  mu.  Value index= 5, #Derivatives= 0\n")
    of.write("  theta.  Value index= 6, #Derivatives= 0\n")

    for i in range(40):
        line = "Node:   %26d" % (int(i)+1)
        print line
        of.write(line + '\n')
        line = ''
        for lambdaCounter in range(4):
            line += "%26.16e" % localParam_0[index[i]*4 + lambdaCounter]
        print line
        of.write(line + '\n')

        line = '%26.16e' % localParam_1[index[i]]
        print line
        of.write(line + '\n')

        line = '%26.16e' % localParam_2[index[i]]
        print line
        of.write(line + '\n')

    of.close()
    f.close()

    
if __name__ == "__main__":
#    if len(sys.argv == 1):
#        print 'Specify model name'
#        sys.exit(2)

    map0 = read_map('CimModelLVPS4x4_0.map')
    map1 = read_map('CimModelLVPS4x4_1.map')
    map2 = read_map('CimModelLVPS4x4_2.map')
#    for list in map2:
#        print list

#    modelName = sys.argv[1]
#    modelName = 'MIDLIFE_01'
#    numFrames = 1
    modelDir = "model"
    files = os.listdir(modelDir)
    print files
    modelFiles = [f for f in files if f.endswith(".model")]

    for f in modelFiles:
        print f
        convert_to_hermite(os.path.join(modelDir,f), map0, map1, map2)
#    for i in range(numFrames):
#        filename = modelName + '_%d.model' % (i + 1)
#        print filename
#        convert_to_hermite(filename, map0, map1, map2)


