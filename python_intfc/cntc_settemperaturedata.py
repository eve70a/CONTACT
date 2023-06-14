
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_settemperaturedata(ire, icp, imeth, params)
#
# set parameters for the temperature calculation for a contact problem
#  imeth          - type of temperature model used (H-digit)
#  nparam         - number of parameters provided
#  params(nparam) - depending on method that is used
#    0: no temperature calculation,     params = []
#    1: keep old parameters,            params = []
#    3: calculate temperature based on new parameters and steady rolling,
#       params = [bktemp1, heatcp1, lambda1, dens1, bktemp2, heatcp2, lambda2, dens2]
#
# dimensions:  bktemp: [C],  heatcp: [J/kg-C],  lambda: [W/length-C],  dens: [kg/length^3]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_settemperaturedata(ire, icp, imeth, params):

    if (not ire):
        ire = 1
    if (not icp):
        icp = -1; # default: W/R contact, all patches
    if (imeth<=0):
        imeth = 0;

    # convert params to NumPy ndarray
    if (not isinstance(params, np.ndarray)):
        if (isinstance(params, list)):
            params = np.array( params, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_settemperaturedata: invalid params provided.')

    nparam = len(params)

    cntc_dll.cntc_settemperaturedata(c_int(ire), c_int(icp), c_int(imeth), c_int(nparam),
                                                                 params.ctypes.data_as(POINTER(c_double)))

# end function cntc_settemperaturedata

#------------------------------------------------------------------------------------------------------------

