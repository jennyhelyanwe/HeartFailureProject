__author__ = 'zwan145'

import os
import numpy


def readOpstre(filename):
    mps = [] # Maximum principal stress at gauss points. 
    fs_2PK = [] # 2PK fibre stress at gauss points. 
    fs_C = [] # Cauchy fibre stress at gauss points. 
    e_idx = [] # Elements indices at gauss points. 
    xi = []  # Xi locations of gauss points. 
    x = []    # Global locations of gauss points. 
    with open(filename, 'r') as file:
        data = file.readline()
        data = file.readline()
        data = file.readline()
        data = file.readline()
        while data != '':
            if data.split()[0] == 'Element':
                e = int(data.split()[1])
                e_idx.append(int(data.split()[1]))
                mps_e = []  # Maximum principal stress per element
                fs_2PK_e = []   # 2PK fibre stress per element
                fs_C_e = []     # Cauchy fibre stress per element
                xi_e = []   # Xi locations of gauss points per element
                x_e = []    # Global coordinates and gauss points per element. 
                for i in range(0,64):
                    x_coord = []
                    temp = file.readline()
                    temp = file.readline()
                    temp = file.readline()
                   
                    # Read global position of gauss point.
                    x_coord.append(float(temp.split()[2]))
                    x_coord.append(float(temp.split()[3]))
                    x_coord.append(float(temp.split()[4]))
                    x_e.append(x_coord)
                    xi_coord = []
                    temp = file.readline()
                    # Read xi location of guass point. 
                    xi_coord.append(float(temp.split()[1]))
                    xi_coord.append(float(temp.split()[2]))
                    xi_coord.append(float(temp.split()[3]))
                    xi_e.append(xi_coord)
                    # Read 2PK fibre stress
                    junk = file.readline()
                    temp = file.readline()
                    fs_2PK_e.append(temp.split()[1])
                    # Read Cauchy fibre stress
                    fs_C_e.append(temp.split()[3])
                    for i in range(0,3):
                        junk = file.readline()
                    # Read maximum principal stress.
                    temp = file.readline()
                    mps_e.append(temp.split()[1])
                    for i in range(0,2):
                        temp = file.readline()
                mps.append(mps_e)
                fs_2PK.append(fs_2PK_e)
                fs_C.append(fs_C_e)
                xi.append(xi_e)
                x.append(x_e)
            data = file.readline()
            data = file.readline()
    return x, xi, mps, fs_2PK, fs_C

def readActiveOpstre(filename):
    mps = []   # Maximum principal stress
    fs_2PK = [] # 2PK fibre stress
    fs_C = []  # Cauchy fibre stress
    e_idx = []  # Element indices
    xi = []  # Xi locations of gauss points. 
    x = []  # Global coordinates of gauss points. 
    with open(filename, 'r') as file:
        data = file.readline()
        data = file.readline()
        data = file.readline()
        data = file.readline()
        while data != '':
            if data.split()[0] == 'Element':
                e_idx.append(int(data.split()[1]))
                mps_e = []
                fs_2PK_e = []
                fs_C_e = []
                xi_e = []
                x_e = []
                for i in range(0,64):
                    x_coord = []
                    temp = file.readline()
                    temp = file.readline()
                    temp = file.readline()
                    # Read global position of gauss point in current element.
                    x_coord.append(float(temp.split()[2]))
                    x_coord.append(float(temp.split()[3]))
                    x_coord.append(float(temp.split()[4]))
                    x_e.append(x_coord)
                    xi_coord = []
                    temp = file.readline()
                    # Read xi position of gauss point in current element. 
                    xi_coord.append(float(temp.split()[1]))
                    xi_coord.append(float(temp.split()[2]))
                    xi_coord.append(float(temp.split()[3]))
                    xi_e.append(xi_coord)
                    # Read 2PK fibre stress
                    junk = file.readline()
                    temp = file.readline()
                    fs_2PK_e.append(temp.split()[1])
                    # Read Cauchy fibre stress
                    fs_C_e.append(temp.split()[3])
                    for i in range(0,3):
                        junk = file.readline()
                    temp= file.readline()
                    # Read maximum principal active stress. 
                    mps_e.append(temp.split()[1])
                    for i in range(0,3):
                        temp = file.readline()
                mps.append(mps_e)
                fs_2PK.append(fs_2PK_e)
                fs_C.append(fs_C_e)
                xi.append(xi_e)
                x.append(x_e)
            data = file.readline()
            data = file.readline()
    return x, xi, mps, fs_2PK, fs_C

def readOpstra(filename):
    exr = []  # Extension ratio.
    e_idx = [] # Elements indices at gauss points.
    xi = []  # Xi locations of gauss points.
    x = []    # Global locations of gauss points.
    with open(filename, 'r') as file:
        for i in range(0,4):
            data = file.readline()
        while data != '':
            if data.split()[0] == 'Element':
                e_idx. append(int(data.split()[1]))
                exr_e = [] # Extension ratio per element
                xi_e = []   # Xi locations of gauss points per element
                x_e = []    # Global coordinates and gauss points per element.
                for i in range(0,64):
                    x_coord = []
                    for i in [1,2,3]:
                        temp = file.readline()
                    #Read global position of gauss point
                    x_coord.append(float(temp.split()[2]))
                    x_coord.append(float(temp.split()[3]))
                    x_coord.append(float(temp.split()[4]))
                    x_e.append(x_coord)
                    xi_coord = []
                    temp = file.readline()
                    # Read xi location of guass point.
                    xi_coord.append(float(temp.split()[1]))
                    xi_coord.append(float(temp.split()[2]))
                    xi_coord.append(float(temp.split()[3]))
                    xi_e.append(xi_coord)
                    junk = file.readline()
                    temp = file.readline()
                    exr_e.append(temp.split()[3])
                    for i in range(0,10):
                        junk = file.readline()
                exr.append(exr_e)
            data = file.readline()
            data = file.readline()
    return x, xi, exr

def readExnode(filename):
    nodes = []
    with open(filename, 'r') as file:
        # Read all junk up till node 9.
        for i in range (0, 151):
            temp = file.readline()
        node_num = int(temp.split()[1])
        # Read node position data.
        while node_num <= 12:
            print node_num
            temp = file.readline()
            nodes_coord_x = []
            for j in range (0, 5):
                nodes_coord_x.append(float(temp.split()[j]))
            temp = file.readline()
            for j in range (0, 3):
                nodes_coord_x.append(float(temp.split()[j]))

            temp = file.readline()
            nodes_coord_y = []
            for j in range (0, 5):
                nodes_coord_y.append(float(temp.split()[j]))
            temp = file.readline()
            for j in range (0, 3):
                nodes_coord_y.append(float(temp.split()[j]))

            temp = file.readline()
            nodes_coord_z = []
            for j in range (0, 5):
                nodes_coord_z.append(float(temp.split()[j]))
            temp = file.readline()
            for j in range (0, 3):
                nodes_coord_z.append(float(temp.split()[j]))

            for j in range(0,7):
                junk = file.readline()
            temp = file.readline()
            node_num = int(temp.split()[1])
            nodes_coord = [nodes_coord_x, nodes_coord_y, nodes_coord_z]
            nodes.append(nodes_coord)
        for i in range(0, 233):
            temp = file.readline()

        node_num = int(temp.split()[1])
        # Read node position data.
        while node_num <= 30:
            print node_num
            temp = file.readline()
            nodes_coord_x = []
            for j in range (0, 5):
                nodes_coord_x.append(float(temp.split()[j]))
            temp = file.readline()
            for j in range (0, 3):
                nodes_coord_x.append(float(temp.split()[j]))

            temp = file.readline()
            nodes_coord_y = []
            for j in range (0, 5):
                nodes_coord_y.append(float(temp.split()[j]))
            temp = file.readline()
            for j in range (0, 3):
                nodes_coord_y.append(float(temp.split()[j]))

            temp = file.readline()
            nodes_coord_z = []
            for j in range (0, 5):
                nodes_coord_z.append(float(temp.split()[j]))
            temp = file.readline()
            for j in range (0, 3):
                nodes_coord_z.append(float(temp.split()[j]))

            for j in range(0,7):
                junk = file.readline()
            temp = file.readline()
            node_num = int(temp.split()[1])

            nodes_coord = [nodes_coord_x, nodes_coord_y, nodes_coord_z]
            nodes.append(nodes_coord)

    return nodes

def readExdata(filename):
    xi = []
    elem = []
    value = []
    with open(filename, 'r') as f:
        for i in range(0, 17):
            junk = f.readline()
        temp = f.readline
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

    return xi, elem, value
