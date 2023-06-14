
#------------------------------------------------------------------------------------------------------------
# function [ pn, px, py ] = cntc_gettractions(ire, icp)
#
# return the tractions for all elements in the potential contact area for a contact problem
#            note the order of the arguments, with pn (z-direction) occurring before px,py.
#
#  mx, my                            - number of elements in potential contact area
#  pn(my,mx), px(my,mx), py(my,mx)   - surface tractions for all elements of contact area, [force/area]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER
from .cntc_getnumelements  import cntc_getnumelements

def cntc_gettractions(ire, icp):

    if (not ire):
        ire = 1
    if (not icp):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_gettractions: not available for icp=%d' % icp)

    mx, my = cntc_getnumelements(ire, icp)

    lenarr = mx * my
    pn = np.zeros(lenarr, dtype=c_double)
    px = np.zeros(lenarr, dtype=c_double)
    py = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_gettractions(c_int(ire), c_int(icp), c_int(lenarr), pn.ctypes.data_as(POINTER(c_double)),
                               px.ctypes.data_as(POINTER(c_double)), py.ctypes.data_as(POINTER(c_double)))

    pn = np.reshape(pn, (my, mx))
    px = np.reshape(px, (my, mx))
    py = np.reshape(py, (my, mx))

    return pn, px, py

# end function cntc_gettractions

#------------------------------------------------------------------------------------------------------------

