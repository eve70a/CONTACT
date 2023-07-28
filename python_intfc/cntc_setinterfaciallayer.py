
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setinterfaciallayer(ire, icp, imeth, params)
#
# set parameters for the interfacial layer for a contact problem
#                                                       *Obsolete, replaced by cntc_setmaterialparameters.
#  imeth          - type of interfacial layer used
#  nparam         - number of parameters provided
#  params(nparam) - depending on method that is used
#
#    0: clean interface, no layer     params = [ ]
#  2,3: Modified FASTSIM algorithm    params = [ k0_mf, alfamf, betamf ]
#    4: elasto-plastic layer          params = [ G3, laythk, tau_c0, k_tau ]
#    
#  dimensions: k0_mf, alfamf, betamf: [-],
#              G3, tau_c0: [force/area],  laythk: [length],  k_tau: [force/area/length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_setinterfaciallayer(ire, icp, imeth, params):
    # default: W/R contact, all patches

    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp = -1
    if (imeth==0):
        params = np.zeros(1, dtype=c_double)
        nparam = 1

    # convert params to NumPy ndarray
    if (not isinstance(params, np.ndarray)):
        if (isinstance(params, list)):
            params = np.array( params, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_setinterfaciallayer: invalid params provided.')

    if (imeth!=0):
        nparam = len(params)

    cntc_dll.cntc_setinterfaciallayer(c_int(ire), c_int(icp), c_int(imeth), c_int(nparam),
                                                                 params.ctypes.data_as(POINTER(c_double)))

# end function cntc_setinterfaciallayer

#------------------------------------------------------------------------------------------------------------

