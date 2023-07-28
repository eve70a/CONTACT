
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setreferencevelocity(ire, icp, veloc)
#
# set the rolling velocity for a contact problem
#
#  veloc          - absolute rolling velocity [veloc]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 1: m=3, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_setreferencevelocity(ire, icp, veloc):

    if (not isinstance(ire, int)):
        ire = 1
    if (not isinstance(icp, int)):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_setreferencevelocity: not available for icp=%d' % icp)

    cntc_dll.cntc_setreferencevelocity(c_int(ire), c_int(icp), c_double(veloc))

# end function cntc_setreferencevelocity

#------------------------------------------------------------------------------------------------------------

