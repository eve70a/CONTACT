
#------------------------------------------------------------------------------------------------------------
# function [ veloc ] = cntc_getreferencevelocity(ire, icp)
#
# get the rolling velocity for a contact problem
#
#  veloc          - absolute rolling velocity [veloc]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_getreferencevelocity(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getreferencevelocity: not available for icp=%d' % icp)

    veloc = c_double()

    cntc_dll.cntc_getreferencevelocity(c_int(ire), c_int(icp), veloc)

    return veloc.value

# end function cntc_getreferencevelocity

#------------------------------------------------------------------------------------------------------------

