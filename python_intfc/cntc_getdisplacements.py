
#------------------------------------------------------------------------------------------------------------
# function [ un, ux, uy ] = cntc_getdisplacements(ire, icp)
#
# return the displ.differences for all elements in the potential contact area for a contact problem
#  mx, my                            - number of elements in potential contact area
#  un(my,mx), ux(my,mx), uy(my,mx)   - displacement difference [length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER
from .cntc_getnumelements  import cntc_getnumelements

def cntc_getdisplacements(ire=1, icp=1):
    print('function called for first time:', __name__)

    if (icp<=0):
        sys.exit('ERROR in cntc_getdisplacements: not available for icp=%d' % icp)

    mx, my = cntc_getnumelements(ire, icp)

    lenarr = mx * my
    un = np.zeros(lenarr, dtype=c_double)
    ux = np.zeros(lenarr, dtype=c_double)
    uy = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getdisplacements(c_int(ire), c_int(icp), c_int(lenarr), un.ctypes.data_as(POINTER(c_double)),
                                    ux.ctypes.data_as(POINTER(c_double)), uy.ctypes.data_as(POINTER(c_double)))

    un = np.reshape(un, (my, mx))
    ux = np.reshape(ux, (my, mx))
    uy = np.reshape(uy, (my, mx))

    return un, ux, uy

# end function cntc_getdisplacements

#------------------------------------------------------------------------------------------------------------

