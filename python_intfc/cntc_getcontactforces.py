
#------------------------------------------------------------------------------------------------------------
# function [ fn, tx, ty, mz ] = cntc_getcontactforces(ire, icp)
#
# return the total forces and torsional moment for a contact problem in contact local coordinates
#
#  fn            - total normal force [force]
#  tx, ty        - total tangential forces [force]
#  mz            - total torsional moment [force.length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_getcontactforces(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getcontactforces: not available for icp=%d' % icp)

    fn = c_double()
    tx = c_double()
    ty = c_double()
    mz = c_double()

    cntc_dll.cntc_getcontactforces(c_int(ire), c_int(icp), fn, tx, ty, mz)

    return fn.value, tx.value, ty.value, mz.value

# end function cntc_getcontactforces

#------------------------------------------------------------------------------------------------------------

