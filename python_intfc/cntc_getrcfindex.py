
#------------------------------------------------------------------------------------------------------------
# function [ rcf ] = cntc_getrcfindex(ire, icp)
#
# return the RCF index for all elements in the potential contact area for a contact problem
#
#  mx, my      - number of elements in potential contact area
#  rcf(my,mx)  - RCF index for all elements of contact area      [-]
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

def cntc_getrcfindex(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getrcfindex: not available for icp=%d' % icp)

    mx, my = cntc_getnumelements(ire, icp)

    lenarr = mx * my
    rcf = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getrcfindex(c_int(ire), c_int(icp), c_int(lenarr), rcf.ctypes.data_as(POINTER(c_double)))

    rcf = np.reshape(rcf, (my, mx))

    return rcf

# end function cntc_getrcfindex

#------------------------------------------------------------------------------------------------------------
