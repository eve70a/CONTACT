
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setextrarigidslip(ire, icp, wx, wy)
#
# set the extra term of the rigid slip for all elements in the potential contact area for a contact problem
#
#  mx, my               - number of elements in potential contact area
#  wx(my,mx), wy(my,mx) - in rolling, T=2,3: extra relative rigid slip  [-]
#                         in shifts,  T=1:   extra rigid shift distance [length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 1: m=3, cp     - require icp>0, default 1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER
from .cntc_getnumelements  import cntc_getnumelements

def cntc_setextrarigidslip(ire, icp, wx, wy):

    if (not ire):
        ire = 1
    if (not icp):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_setextrarigidslip: not available for icp=%d' % icp)

    mx, my = cntc_getnumelements(ire, icp)

    lenarr = mx * my

    # reshape to 1D array, with ix running fastest(?)
    if (np.ndim(wx)>1):
        if (not (np.shape(wx) == (my, mx))):
            sys.exit('ERROR: wx should be 1D or (my,mx) array')
        wx = np.reshape(wx, (lenarr))
    if (np.ndim(wy)>1):
        if (not np.shape(wy) == (my, mx)):
            sys.exit('ERROR: wy should be 1D or (my,mx) array')
        wy = np.reshape(wy, (lenarr))

    cntc_dll.cntc_setextrarigidslip(c_int(ire), c_int(icp), c_int(lenarr), 
                                   wx.ctypes.data_as(POINTER(c_double)), wy.ctypes.data_as(POINTER(c_double)))

# end function cntc_setextrarigidslip

#------------------------------------------------------------------------------------------------------------

