
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setwheelsetposition(ire, ewheel, params) 
#
# set the wheelset position state data for a wheel-rail contact problem
#  ewheel         - type of position specification (E-digit)
#  nparam         - number of parameters provided
#  params(nparam) - depending on method that is used
#
#  E=0  : keep wheelset position from previous specification, ignore params provided
#  E=1-8: new wheelset position   params = [s, y, z, roll, yaw, pitch]
#
#  dimensions:   s_ws, y_ws, z_ws [length],       roll, yaw, pitch [angle]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_setwheelsetposition(ire, ewheel, params):

    if (not isinstance(ire, int)):
        ire = 1

    # convert params to NumPy ndarray
    if (not isinstance(params, np.ndarray)):
        if (isinstance(params, list)):
            params = np.array( params, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_setwheelsetposition: invalid params provided.')

    nparam = len(params)

    cntc_dll.cntc_setwheelsetposition(c_int(ire), c_int(ewheel), c_int(nparam),
                                                                 params.ctypes.data_as(POINTER(c_double)))

# end function cntc_setwheelsetposition

#------------------------------------------------------------------------------------------------------------

