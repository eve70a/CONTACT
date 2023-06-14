
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_settrackdimensions(ire, ztrack, params) 
#
# set the track or roller-rig description for a wheel-rail contact problem
#
#   !!!  This is has become an alias for cntc_settrackdimensions_new.
#        This used to link to cntc_settrackdimensions_old in the previous version. !!!
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_settrackdimensions(ire, ztrack, params):

    if (not ire):
        ire = 1

    # convert params to NumPy ndarray
    if (not isinstance(params, np.ndarray)):
        if (isinstance(params, list)):
            params = np.array( params, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_settrackdimensions: invalid params provided.')

    nparam = len(params)

    cntc_dll.cntc_settrackdimensions(c_int(ire), c_int(ztrack), c_int(nparam), 
                                                                 params.ctypes.data_as(POINTER(c_double)))

# end function cntc_settrackdimensions

#------------------------------------------------------------------------------------------------------------

