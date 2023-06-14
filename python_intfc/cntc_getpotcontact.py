
#------------------------------------------------------------------------------------------------------------
# function [ mx, my, xc1, yc1, dx, dy ] = cntc_getpotcontact(ire, icp)
#
# get the parameters of the potential contact area for a contact problem
#    3: first center + grid sizes,        params = [ mx, my, xc1, yc1, dx, dy ]
#
#  mx, my        - number of elements in x- and y-directions [-]
#  xc1, yc1      - position of first element center [length]
#  dx, dy        - grid discretization step sizes [length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getpotcontact(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getpotcontact: not available for icp=%d' % icp)

    lenarr = 6
    values = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getpotcontact(c_int(ire), c_int(icp), c_int(lenarr), 
                                                                 values.ctypes.data_as(POINTER(c_double)))

    mx  = int(values[0])
    my  = int(values[1])
    xc1 = values[2]
    yc1 = values[3]
    dx  = values[4]
    dy  = values[5]

    return mx, my, xc1, yc1, dx, dy

# end function cntc_getpotcontact

#------------------------------------------------------------------------------------------------------------
