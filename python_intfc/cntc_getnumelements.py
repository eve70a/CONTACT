
#------------------------------------------------------------------------------------------------------------
# function [ mx, my ] = cntc_getnumelements(ire, icp)
#
# return the number of elements in the potential contact area used for a contact problem,
#            length of tractions arrays
#
#  mx, my        - number of discretization elements in long/lat dirs
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int

def cntc_getnumelements(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getnumelements: not available for icp=%d' % icp)

    mx = c_int()
    my = c_int()

    cntc_dll.cntc_getnumelements(c_int(ire), c_int(icp), mx, my)

    return mx.value, my.value

# end function cntc_getnumelements

#------------------------------------------------------------------------------------------------------------

