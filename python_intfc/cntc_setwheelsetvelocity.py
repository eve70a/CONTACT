
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setwheelsetvelocity(ire, ewheel, params) 
#
# continue the wheelset description for a wheel-rail contact problem
#  ewheel         - type of velocity specification (E-digit)
#  nparam         - number of parameters provided
#  params(nparam) - depending on method that is used
#
#  E=0-1: keep velocity settings from previous specification, ignore params provided
#  E=2-8: new wheelset velocity   params = [vs, vy, vz, vroll, vyaw, vpitch]
#
#  dimensions:  vs_ws, vy_ws, vz_ws [veloc],     vroll, vyaw, vpitch [ang.veloc]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_setwheelsetvelocity(ire, ewheel, params):

    if (not isinstance(ire, int)):
        ire = 1

    # convert params to NumPy ndarray
    if (not isinstance(params, np.ndarray)):
        if (isinstance(params, list)):
            params = np.array( params, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_setwheelsetvelocity: invalid params provided.')

    nparam = len(params)

    cntc_dll.cntc_setwheelsetvelocity(c_int(ire), c_int(ewheel), c_int(nparam),
                                                                 params.ctypes.data_as(POINTER(c_double)))

# end function cntc_setwheelsetvelocity

#------------------------------------------------------------------------------------------------------------

