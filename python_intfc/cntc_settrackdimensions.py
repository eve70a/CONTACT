
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_settrackdimensions(ire, ztrack, params) 
#
# set the track or roller-rig description for a wheel-rail contact problem
#  ztrack    - control digit ZTRACK
#  params    - depending on method that is used
#
#    0: maintain track dimensions     params = [ ]
#    1: new design track dimensions   params = [gaught, gaugsq, gaugwd, cant, nomrad],   if gaught >  0,
#                                         or   [gaught, raily0, railz0, cant, nomrad],   if gaught <= 0.
#    2: new track deviations          params = [dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
#    3: new dimensions & track deviations for current side of the track
#                                     params = params(1:5) cf. Z=1 followed by params(6:11) cf. Z=2;
#                                              additionally, [kyrail, fyrail, kzrail, fzrail] when F=3.
#
# dimensions: gaught, gaugwd, raily0, railz0, nomrad, dyrail, dzrail [length],    cant, drollr [angle]
#                                                     vyrail, vzrail [veloc],           vrollr [ang.veloc]
#                                                kyrail, kzrail [force/length], fyrail, fzrail [force]
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

    if (not isinstance(ire, int)):
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

