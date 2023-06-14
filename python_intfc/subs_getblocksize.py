
#------------------------------------------------------------------------------------------------------------
# function [ nx, ny, nz ] = subs_getblocksize(ire, icp, iblk)
#
# get the number of points in a block used for subsurface stress calculation
#
#  nx, ny, nz    - number of points used in x-, y- and z-directions
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def subs_getblocksize(ire=1, icp=1, iblk=1):

    if (icp<=0):
        sys.exit('ERROR in subs_getresults: not available for icp=%d' % icp)

    nx = c_int()
    ny = c_int()
    nz = c_int()

    cntc_dll.subs_getblocksize(c_int(ire), c_int(icp), c_int(iblk), nx, ny, nz)

    return nx.value, ny.value, nz.value

# end function subs_getblocksize

#------------------------------------------------------------------------------------------------------------
