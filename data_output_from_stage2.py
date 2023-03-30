import yt                                                                                                       
import numpy as np
import h5py
import time
from enum import Enum
from global_variables import *


def SetOutflowBoundary(BaryonicField, N_ghost=3):
    # obtain the array shape of the datadump. Ex: [256, 256, 256]
    SizeofDD         = np.array(BaryonicField.shape)
    print("Original shape = {}".format(SizeofDD))
    # add the number of ghost zones
    SizeofDD        += 2*N_ghost
    print("After adding ghost zones = {}".format(SizeofDD))
    # create new data array with size including ghost zones. [262, 262, 262]
    ArrayDD          = np.zeros(SizeofDD)
    # place the data dump into the new array at the corresponding location. [3~259, 3~259, 3~259]
    ArrayDD[N_ghost:-N_ghost, N_ghost:-N_ghost, N_ghost:-N_ghost] = BaryonicField
    # set the rest zones according to the outflow boundary condition q(-x)=q(0).
    for i in range(N_ghost-1, -1, -1): 
        for j in range(SizeofDD[0]):
            for k in range(SizeofDD[0]):
                ArrayDD[i,j,k] = ArrayDD[i+1,j,k] # frount x axis: [2~0, :, :] 
                ArrayDD[j,k,i] = ArrayDD[j,k,i+1] # lower z axis: [:, :, 2~0]
                ArrayDD[k,i,j] = ArrayDD[k,i+1,j] # left y axis: [:, 2~0, :]

    for i in range(SizeofDD[0]-N_ghost, SizeofDD[0]):
        for j in range(SizeofDD[0]):
            for k in range(SizeofDD[0]):
                ArrayDD[i,j,k] = ArrayDD[i-1,j,k] # x axis: [260~262, :, :]
                ArrayDD[j,k,i] = ArrayDD[j,k,i-1] # z axis: [:, :, 260~262]
                ArrayDD[k,i,j] = ArrayDD[k,i-1,j] # y axis: [:, 260~262, :]

    return ArrayDD


def DataOutputWithNewBoundary(file_dir, output_dir, dumpnum, SetTheBoundary, N_ghost=3):
    dataset           = yt.load("{}/DD{:0>4d}/DD{:0>4d}".format(file_dir, dumpnum, dumpnum))
    # access specific data 
    all_data_level_0  = dataset.covering_grid(level=0, left_edge=[0.0,0.0,0.0], dims=dataset.domain_dimensions)

    # Now we open our output file using h5py
    f                 = h5py.File(output_dir, mode="w")
    # List of target fields
    TargetFields      = ["Density","GasEnergy","TotalEnergy","x-velocity","y-velocity","z-velocity",
                         "Electron_Density","H2II_Density","H2I_Density","HII_Density","HI_Density","HM_Density",
                         "HeIII_Density","HeII_Density","HeI_Density"]
     
    # Rewrite the field data with outflow Boundary
    for field in TargetFields:
        print("Setting new boundary for {}.".format(field))
        DataForStage3 = SetTheBoundary(all_data_level_0[field], N_ghost)
        f.create_dataset("/"+field, data=DataForStage3)
        print("Field {0} has been stored with new outflow boundary condition. (shape={1})\n".format(field, DataForStage3.shape))
    
    # We close the file
    f.close()


