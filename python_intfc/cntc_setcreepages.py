
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setcreepages(ire, icp, vx, vy, phi)
#
# set the kinematic constants (creepages) for a contact problem
# note: vx is ignored when the F-digit is 1 or 2, vy is ignored when F=1.
#
#  vx, vy, phi    - in rolling, T=2,3: long/lat/spin creepages [-, -, angle/length]
#                   in shifts,  T=1:   long/lat/spin shift     [length, length, angle]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 1: m=3, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_setcreepages(ire, icp, vx, vy, phi):

    if (not ire):
        ire = 1
    if (not icp):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_setcreepages: not available for icp=%d' % icp)

    cntc_dll.cntc_setcreepages(c_int(ire), c_int(icp), c_double(vx), c_double(vy), c_double(phi))

# end function cntc_setcreepages

#------------------------------------------------------------------------------------------------------------

