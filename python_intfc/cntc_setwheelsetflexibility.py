
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setwheelsetflexibility(ire, ewheel, params) 
#
# set the description of wheelset flexibilities for a wheel-rail contact problem
#    ewheel         - type of wheelset flexibilities specification (E-digit)
#    nparam         - number of parameters provided
#    params(nparam) - depending on method that is used
#
#    E=0  : keep wheelset flexibility from previous specification, ignore params provided
#    E=1-4: no wheelset flexibility               params = []
#    E=5,7: new wheelset flexibility parameters   params = [dxwhl, dywhl, dzwhl, drollw, dyaww, dpitchw,
#                                                           vxwhl, vywhl, vzwhl, vrollw, vyaww, vpitchw],
#    E=6,8: same as 5,7, with separate values for both sides of the wheelset
#
#    dimensions:   dxwhl, dywhl, dzwhl [length],  drollw, dyaww, dpitchw [angle]
#                  vxwhl, vywhl, vzwhl [veloc],   vrollw, vyaww, vpitchw [ang.veloc]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_setwheelsetflexibility(ire, ewheel, params):

    if (not isinstance(ire, int)):
        ire = 1
    if (not isinstance(ewheel, int)):
        ewheel = 0
    if (ewheel<=4):
        params = np.zeros(1, dtype=c_double)

    # convert params to NumPy ndarray
    if (not isinstance(params, np.ndarray)):
        if (isinstance(params, list)):
            params = np.array( params, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_setwheelsetflexibility: invalid params provided.')

    nparam = len(params)

    cntc_dll.cntc_setwheelsetflexibility(c_int(ire), c_int(ewheel), c_int(nparam),
                                                                 params.ctypes.data_as(POINTER(c_double)))

# end function cntc_setwheelsetflexibility

#------------------------------------------------------------------------------------------------------------

