
#------------------------------------------------------------------------------------------------------------
# function [ tcpu, twall ] = cntc_getcalculationtime(ire, icp)
#
# return accumulated cpu-time and wall-clock-time used since last timer reset for a contact problem
#  tcpu, twall   - cpu- and wall-clock times used [time]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, byref, POINTER

def cntc_getcalculationtime(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getcalculationtime: not available for icp=%d' % icp)

    tcpu  = c_double(-1)
    twall = c_double(-9)

    cntc_dll.cntc_getcalculationtime(c_int(ire), c_int(icp), tcpu, twall)

    return tcpu.value, twall.value

# end function cntc_getcalculationtime

#------------------------------------------------------------------------------------------------------------

