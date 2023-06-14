
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setnormalforce(ire, icp, fn)
#
# set the total normal force of the bodies as a whole for a contact problem
# Note: this function sets control digit N = 1
#
#  fn             - total normal force between the two bodies [force]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 1: m=3, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_setnormalforce(ire, icp, fn):

    if (not ire):
        ire = 1
    if (not icp):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_setnormalforce: not available for icp=%d' % icp)

    cntc_dll.cntc_setnormalforce(c_int(ire), c_int(icp), c_double(fn))

# end function cntc_setnormalforce

#------------------------------------------------------------------------------------------------------------

