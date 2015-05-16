## This script aligns the ED epicardial and endocardial surface data 
# with the current model ED prediction. This allows the surface data to be
# projected to the ED model and give error values which are only sensitive 
# to material properties changes and relatively insensitive to erroneous 
# basal displacement boundary conditions. 

# This script should read in the current ED geometry in an ipnode file, and also the 
# ED surface data in ipdata file and calculate the transformation matrix required for
# alignment. The actual rotation must then be performed inside the CMISS command file 
# OptimiseActive_Holzapfel_M.com before the projection is performed. 

# Generate surface data for current ED model

import os
import scipy
import scipy.constants
import numpy
import math
import fitting_SurfaceData
import re
import string

from scipy import array,concatenate
from scipy import spatial
from scipy import delete
from fitting_SurfaceData import readSurfaceData, readIpdata, constructTransformationMatrix, inverseTransformRigid3D, transformScale3D
from fitting_SurfaceData import transformRigid3DFinal, fitDataRigidScaleNoCorr, fitDataRigidScaleNoCorr_TwoSurfaces, fitDataRigidAnisotropicScaleNoCorr
from fitting_SurfaceData import fitDataRigidScaleNoCorr_ModelTree, fitDataRigidAnisotropicScaleNoCorr_ModelTree, fitDataRigidAnisotropicScaleNoCorr_ModelTree_TwoSurfaces
from fitting_SurfaceData import fitDataAnisotropicScaleRigidNoCorr_ModelTree_TwoSurfaces
from fitting_SurfaceData import writeIpdata, writeTransformation, writeTransformation_Scaling, writeTransformation_Rotation


# Read in current ED model surface data
filename_Epi = 'output_debug/ES_model_surface_epi.ipdata'
filename_Endo = 'output_debug/ES_model_surface_endo.ipdata'

modelEpi = readIpdata(filename_Epi)
nd_modelEpi = len(modelEpi)
modelEndo = readIpdata(filename_Endo)
nd_modelEndo = len(modelEndo)

model = concatenate((modelEpi, modelEndo),0)
nd_model = len(model)

print 'Total number of data from the model is ', nd_model

## Read in ED experimental surface data
SDEpi = readIpdata('Surface_Points_Epi_ES.ipdata')
SDEndo = readIpdata('Surface_Points_Endo_ES.ipdata')
SD = concatenate((SDEpi, SDEndo),0)
nd_SD = len(SD)
print 'Total number of experimental data is ', nd_SD

[transVector, SDEpi_transformed, SDEndo_transformed] = fitDataAnisotropicScaleRigidNoCorr_ModelTree_TwoSurfaces(SDEpi, SDEndo, modelEpi, modelEndo, 'output_debug/TransformationMatrix.TRN', xtol=1e-5, maxfev=0, t0=None)

print 'Transformation vector is'
print transVector
print ''


