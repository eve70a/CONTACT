
#------------------------------------------------------------------------------------------------------------
# function [ pen ] = cntc_getpenetration(ire, icp)
#
# return the penetration (approach) for a contact problem
#
#  pen           - penetration [length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_getpenetration(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getpenetration: not available for icp=%d' % icp)

    pen = c_double()

    cntc_dll.cntc_getpenetration(c_int(ire), c_int(icp), pen)

    return pen.value

# end function cntc_getpenetration

#------------------------------------------------------------------------------------------------------------

