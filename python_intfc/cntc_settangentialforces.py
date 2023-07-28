
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_settangentialforces(ire, icp, fx, fy)
#
# set the total tangential forces a contact problem
# note: fx is ignored when the F-digit is 0, fy is ignored when F=0 or 1.
#
#  fx, fy       - total tangential forces relative to fstat*fn [-]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 1: m=3, cp     - require icp>0, default 1

from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_settangentialforces(ire, icp, fx, fy):

    if (not isinstance(ire, int)):
        ire = 1
    if (not isinstance(icp, int)):
        icp = 1

    cntc_dll.cntc_settangentialforces(c_int(ire), c_int(icp), c_double(fx), c_double(fy))

# end function cntc_settangentialforces

#------------------------------------------------------------------------------------------------------------

