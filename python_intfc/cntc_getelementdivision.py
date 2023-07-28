
#------------------------------------------------------------------------------------------------------------
# function [ eldiv ] = cntc_getelementdivision(ire, icp)
#
# return flags for all elements in the potential contact area for a contact problem
# indicating whether the element is in Exerior (0), Adhesion (1), Slip (2) or Plasticity (3).
#
#  mx, my       - number of elements in potential contact area
#  eldiv(my,mx) - element division of contact area, 0=E, 1=H, 2=S, 3=P.
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

def cntc_getelementdivision(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getelementdivision: not available for icp=%d' % icp)

    mx, my = cntc_getnumelements(ire, icp)

    lenarr = mx * my
    eldiv  = np.zeros(lenarr, dtype=c_int)

    cntc_dll.cntc_getelementdivision(c_int(ire), c_int(icp), c_int(lenarr), 
                                                                       eldiv.ctypes.data_as(POINTER(c_int)))

    # TODO: convert to double for consistency with loadcase.m ?
    eldiv = np.reshape(eldiv, (my, mx))

    return eldiv

# end function cntc_getelementdivision

#------------------------------------------------------------------------------------------------------------

