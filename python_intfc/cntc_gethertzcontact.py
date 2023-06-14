
#------------------------------------------------------------------------------------------------------------
# function [ rvalues ] = cntc_gethertzcontact(ire, icp)
#
# get the parameters from a Hertzian contact problem
#
# The following values are returned in rvalues:
#   1 - A1       - curvature in rolling direction [1/length]
#   2 - B1       - curvature in lateral direction [1/length]
#   3 - AA       - semi-axis in rolling direction [length]
#   4 - BB       - semi-axis in lateral direction [length]
#   5 - RHO      - effective radius of curvature, 2 / (A1 + B1) [length]
#   6 - CP       - effective semi-axis, sqrt(AA * BB) [length]
#   7 - SCALE    - potential contact scale factor [-]
#   8 - BNEG     - semi-axis of negative half-ellipse in lateral direction [length]
#   9 - BPOS     - semi-axis of positive half-ellipse in lateral direction [length]
#  10 - AOB      - ellipticity AA/BB [-]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 1: m=3, cp     - require icp>0, default 1

import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_gethertzcontact(ire=1, icp=-1):

    lenarr = 10
    values = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_gethertzcontact(c_int(ire), c_int(icp), c_int(lenarr), 
                                                                 values.ctypes.data_as(POINTER(c_double)))

    return values

# end function cntc_getglobalforces

#------------------------------------------------------------------------------------------------------------
