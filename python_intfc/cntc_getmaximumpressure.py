
#------------------------------------------------------------------------------------------------------------
# function [ pmax ] = cntc_getmaximumpressure(ire, icp)
#
# return the maximum normal pressure in a contact problem
#
#  pnmax         - maximum pressure in contact patch [force/area]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_getmaximumpressure(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getmaximumpressure: not available for icp=%d' % icp)

    pmax = c_double()

    cntc_dll.cntc_getmaximumpressure(c_int(ire), c_int(icp), pmax)

    return pmax.value

# end function cntc_getmaximumpressure

#------------------------------------------------------------------------------------------------------------

