
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setpotcontact(ire, icp, ipotcn, params)
#
# set the parameters of the potential contact area for a contact problem
#  ipotcn         - type of specification for the potential contact area
#  nparam         - number of parameters provided
#  params(nparam) - depending on method that is used
#
#  for w/r contact, icp = -1:
#    0: w/r contact, fixed grid sizes,     params = [ dx, ds, n.a. ]
#   -1: w/r contact, fixed grid sizes,  params = [ dx, ds, a_sep, d_sep, d_comb, [d_turn] ]
#
#  for generic contact, icp > 0:
#    1: lower-left + grid sizes,        params = [ mx, my, xl , yl , dx , dy  ]
#    2: lower-left + upper right,       params = [ mx, my, xl , yl , xh , yh  ]
#    3: 1st center + grid sizes,        params = [ mx, my, xc1, yc1, dx , dy  ]
#    4: 1st center + last center,       params = [ mx, my, xc1, yc1, xcm, ycm ]
#
#  dimensions: mx, my [-],   a_sep [angle],
#              dx, ds, dy, d_sep, d_comb, d_turn, xl, yl, xh, yh, xc1, yc1, xcm, ycm [length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_setpotcontact(ire, icp, ipotcn, params):
    # default: W/R contact, all patches

    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp = -1
    # convert params to NumPy ndarray
    if (not isinstance(params, np.ndarray)):
        if (isinstance(params, list)):
            params = np.array( params, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_setpotcontact: invalid params provided.')

    nparam = len(params)

    cntc_dll.cntc_setpotcontact(c_int(ire), c_int(icp), c_int(ipotcn), c_int(nparam),
                                                                 params.ctypes.data_as(POINTER(c_double)))

# end function cntc_setpotcontact

#------------------------------------------------------------------------------------------------------------

