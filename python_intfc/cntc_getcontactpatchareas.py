
#------------------------------------------------------------------------------------------------------------
# function [ carea, harea, sarea ] = cntc_getcontactpatchareas(ire, icp)
#
# return the area of contact for a contact problem
#  carea   - area of contact patch [area]
#  harea   - area of adhesion area [area]
#  sarea   - area of slip area [area]
#  parea   - area of plasticity area [area]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_getcontactpatchareas(ire=1, icp=1):
    print('function called for first time:', __name__)

    if (icp<=0):
        sys.exit('ERROR in cntc_getcontactpatchareas: not available for icp=%d' % icp)

    carea  = c_double()
    harea  = c_double()
    sarea  = c_double()

    cntc_dll.cntc_getcontactpatchareas(c_int(ire), c_int(icp), carea, harea, sarea)

    parea  = carea.value - harea.value - sarea.value

    return carea.value, harea.value, sarea.value, parea

# end function cntc_getcontactpatchareas

#------------------------------------------------------------------------------------------------------------

