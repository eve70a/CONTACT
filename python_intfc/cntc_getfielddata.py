#------------------------------------------------------------------------------------------------------------
# function [ fld ] = cntc_getfielddata(ire, icp, ifld)
#
# return the values of field 'ifld' for all elements in the potential contact area for a contact problem
#
#  mx, my        - number of elements in potential contact area
#  fld(my,mx)    - output value of field 'ifld' for all elements of contact area
#
#                  ifld  =  CNTC_fld_h      =    1 ! for retrieving array h      [length]
#                           CNTC_fld_mu     =    2 ! for retrieving array mu     [-]
#                           CNTC_fld_px     =    3 ! for retrieving array px     [force/area]
#                           CNTC_fld_py     =    4 ! for retrieving array py     [force/area]
#                           CNTC_fld_pn     =    5 ! for retrieving array pn     [force/area]
#                           CNTC_fld_ux     =    7 ! for retrieving array ux     [length]
#                           CNTC_fld_uy     =    8 ! for retrieving array uy     [length]
#                           CNTC_fld_un     =    9 ! for retrieving array un     [length]
#                           CNTC_fld_taucrt =   11 ! for retrieving array taucrt [force/area]
#                           CNTC_fld_uplsx  =   12 ! for retrieving array uplsx  [length]
#                           CNTC_fld_uplsy  =   13 ! for retrieving array uplsy  [length]
#                           CNTC_fld_sx     =   15 ! for retrieving array sx     [-]
#                           CNTC_fld_sy     =   16 ! for retrieving array sy     [-]
#                           CNTC_fld_temp1  =   20 ! for retrieving array temp1  [C]
#                           CNTC_fld_temp2  =   21 ! for retrieving array temp2  [C]
#                           CNTC_fld_wx     =   22 ! for retrieving array wx     [-]
#                           CNTC_fld_wy     =   23 ! for retrieving array wy     [-]
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

def cntc_getfielddata(ire, icp, ifld):

    if (not isinstance(ire, int)):
        ire = 1
    if (not isinstance(icp, int)):
        icp = 1
    if (not isinstance(ifld, int)):
        ifld = CNTC.fld_pn
    if (icp<=0):
        sys.exit('ERROR in cntc_getfielddata: not available for icp=%d' % icp)

    mx, my = cntc_getnumelements(ire, icp)

    lenarr = mx * my
    val = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getfielddata(c_int(ire), c_int(icp), c_int(ifld), c_int(lenarr), 
                               val.ctypes.data_as(POINTER(c_double)))

    val = np.reshape(val, (my, mx))

    return val

# end function cntc_getfielddata

#------------------------------------------------------------------------------------------------------------

