
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_settrackdimensions(ire, ztrack, params) 
#
# set the track or roller-rig description for a wheel-rail contact problem
#  ztrack - control digit ZTRACK      params - depending on method that is used
#    0: maintain track dimensions     params = [ ]
#    1: new design track dimensions   params = [gaught,  dummy, gaugwd, cant, nomrad, curv],  if gaught >  0,
#                                         or   [gaught, raily0, railz0, cant, nomrad, curv],  if gaught <= 0.
#    2: new track deviations          params = [dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
#    3: new dimensions & track deviations for current side of the track
#                                     params = params(1:6) cf. Z=1 followed by params(7:12) cf. Z=2;
#  ztrack >= 30 is used to configure the massless rail model, F=3
#   32/33: "F=3, Z=2/3":              params = [kyrail, fyrail, kzrail, fzrail, {dystep0} ]
#
# dimensions: gaught, gaugwd, raily0, railz0, nomrad  [length],  cant  [angle],  curv [1/length]
#             dyrail, dzrail  [length],  drollr  [angle],  vyrail, vzrail  [veloc], vrollr  [ang.veloc]
#             kyrail, kzrail  [force/length],  fyrail, fzrail  [force],  dystep0  [length]
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
