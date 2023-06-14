
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_settimestep(ire, icp, dt)
#
# set the time step size dt for a contact problem, in particular for T = 0 or 1 (shifts)
#
#  dt          - time step size [time]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 1: m=3, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_settimestep(ire, icp, dt):

    if (not ire):
        ire = 1
    if (not icp):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_settimestep: not available for icp=%d' % icp)

    cntc_dll.cntc_settimestep(c_int(ire), c_int(icp), c_double(dt))

# end function cntc_settimestep

#------------------------------------------------------------------------------------------------------------

