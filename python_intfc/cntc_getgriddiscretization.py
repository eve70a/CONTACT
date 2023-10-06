
#------------------------------------------------------------------------------------------------------------
# function [ dx, dy ] = cntc_getgriddiscretization(ire, icp)
#
# get the grid discretization step sizes dx,dy for a contact problem
#
#  dx, dy        - grid discretization step sizes [length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_getgriddiscretization(ire=1, icp=1):
    print('function called for first time:', __name__)

    if (icp<=0):
       sys.exit('ERROR in cntc_getgriddiscretization: not available for icp=%d' % icp)

    dx = c_double()
    dy = c_double()

    cntc_dll.cntc_getgriddiscretization(c_int(ire), c_int(icp), dx, dy)

    return dx.value, dy.value

# end function cntc_getgriddiscretization

#------------------------------------------------------------------------------------------------------------

