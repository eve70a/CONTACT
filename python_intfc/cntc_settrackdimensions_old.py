
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_settrackdimensions_old(ire, ztrack, params) 
#
# set the track or roller-rig description for a wheel-rail contact problem
#  ztrack         - control digit ZTRACK
#  nparam         - number of parameters provided
#  params(nparam) - depending on method that is used
#
#    1: new design track dimensions   params = [gaugwd, gaught, cant, nomrad ]
#    2: new track deviations          params = [dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
#    3: new dimensions & deviations for current side of the track
#           params = [gaugwd, gaught, cant, nomrad, dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
#
# dimensions: gaugwd, gaught, nomrad, dyrail, dzrail [length],  cant, drollr [angle],
#                                     vyrail, vzrail [veloc],         vrollr [ang.veloc]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_settrackdimensions_old(ire, ztrack, params):

    if (not isinstance(ire, int)):
        ire = 1

    # convert params to NumPy ndarray
    if (not isinstance(params, np.ndarray)):
        if (isinstance(params, list)):
            params = np.array( params, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_settrackdimensions_old: invalid params provided.')

    nparam = len(params)

    cntc_dll.cntc_settrackdimensions_old(c_int(ire), c_int(ztrack), c_int(nparam), 
                                                                 params.ctypes.data_as(POINTER(c_double)))

# end function cntc_settrackdimensions_old

#------------------------------------------------------------------------------------------------------------

