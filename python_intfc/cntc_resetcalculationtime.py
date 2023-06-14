
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_resetcalculationtime(ire, icp)
#
# reset the accumulated cpu-time and wall-clock-time used for a contact problem
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int

def cntc_resetcalculationtime(ire=1, icp=1):

    if (icp<=0):
       sys.exit('ERROR in cntc_resetcalculationtime: not available for icp=%d' % icp)

    cntc_dll.cntc_resetcalculationtime(c_int(ire), c_int(icp))

# end function cntc_resetcalculationtime

#------------------------------------------------------------------------------------------------------------

