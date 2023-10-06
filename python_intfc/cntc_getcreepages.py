
#------------------------------------------------------------------------------------------------------------
# function [ vx, vy, phi ] = cntc_getcreepages(ire, icp)
#
# get the kinematic constants (creepages) for a contact problem
#
#  vx, vy, phi    - in rolling, T=2,3: long/lat/spin creepages [-, -, angle/length]
#                   in shifts,  T=1:   long/lat/spin shift [length, length, angle]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getcreepages(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getcreepages: not available for icp=%d' % icp)

    vx  = c_double()
    vy  = c_double()
    phi = c_double()

    cntc_dll.cntc_getcreepages(c_int(ire), c_int(icp), vx, vy, phi)

    return vx.value, vy.value, phi.value

# end function cntc_getcreepages

#------------------------------------------------------------------------------------------------------------

