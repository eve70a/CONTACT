
#------------------------------------------------------------------------------------------------------------
# function [ t1max, t2max ] = cntc_getmaximumtemperature(ire, icp)
#
# return the maximum contact temperature in a contact problem
#
#  t1max, t2max  - maximum surface temperatures in bodies 1 and 2 in contact patch [C]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_getmaximumtemperature(ire=1, icp=1):
    print('function called for first time:', __name__)

    if (icp<=0):
        sys.exit('ERROR in cntc_getmaximumtemperature: not available for icp=%d' % icp)

    t1max = c_double()
    t2max = c_double()

    cntc_dll.cntc_getmaximumtemperature(c_int(ire), c_int(icp), t1max, t2max)

    return t1max.value, t2max.value

# end function cntc_getmaximumtemperature

#------------------------------------------------------------------------------------------------------------

