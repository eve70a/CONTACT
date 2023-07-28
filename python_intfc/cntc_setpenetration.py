
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setpenetration(ire, icp, pen)
#
# set the approach (penetration) of the bodies as a whole for a contact problem
# Note: this function sets control digit N = 0
#
#  pen            - penetration/approach of the two bodies [length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 1: m=3, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_setpenetration(ire, icp, pen):

    if (not isinstance(ire, int)):
        ire = 1
    if (not isinstance(icp, int)):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_setpenetration: not available for icp=%d' % icp)

    cntc_dll.cntc_setpenetration(c_int(ire), c_int(icp), c_double(pen))

# end function cntc_setpenetration

#------------------------------------------------------------------------------------------------------------

