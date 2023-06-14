
#------------------------------------------------------------------------------------------------------------
# function [ sx, sy ] = cntc_getmicroslip(ire, icp)
#
# return the relative micro-slip velocity for all elements in the potential contact area for
#            a contact problem
#
#  mx, my               - number of elements in potential contact area
#  sx(my,mx), sy(my,mx) - in rolling, T=2,3: relative micro-slip velocity [-]
#                         in shifts,  T=1:   shift distance [length]
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

def cntc_getmicroslip(ire, icp, mx, my):

    if (not ire):
        ire = 1
    if (not icp):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_getmicroslip: not available for icp=%d' % icp)

    mx, my = cntc_getnumelements(ire, icp)

    lenarr = mx * my
    sx = np.zeros(lenarr, dtype=c_double)
    sy = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getmicroslip(c_int(ire), c_int(icp), c_int(lenarr), 
                              sx.ctypes.data_as(POINTER(c_double)), sy.ctypes.data_as(POINTER(c_double)))

    sx = np.reshape(sx, (my, mx))
    sy = np.reshape(sy, (my, mx))

    return sx, sy

# end function cntc_getmicroslip

#------------------------------------------------------------------------------------------------------------

