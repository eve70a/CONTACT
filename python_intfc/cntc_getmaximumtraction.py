
#------------------------------------------------------------------------------------------------------------
# function [ ptmax ] = cntc_getmaximumtraction(ire, icp)
#
# return the maximum tangential traction in a contact problem
#
#  ptmax         - maximum traction |pt| in contact patch [force/area]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_getmaximumtraction(ire=1, icp=1):
    print('function called for first time:', __name__)

    if (icp<=0):
        sys.exit('ERROR in cntc_getmaximumtraction: not available for icp=%d' % icp)

    ptmax = c_double()

    cntc_dll.cntc_getmaximumtraction(c_int(ire), c_int(icp), ptmax)

    return ptmax.value

# end function cntc_getmaximumtraction

#------------------------------------------------------------------------------------------------------------

