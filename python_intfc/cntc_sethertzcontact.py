
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_sethertzcontact(ire, icp, ipotcn, params)
#
# set the parameters for a Hertzian contact problem
#  ipotcn         - type of specification of the Hertzian geometry
#  nparam         - number of parameters provided
#  params(nparam) - depending on method that is used
#
#   -6: SDEC approach, union of two half ellipses                  params = [ mx, my, aa , bneg, bpos, scale ]
#   -5: Hertzian rectangular contact, half-sizes prescribed,       params = [ mx, my, aa , bb , scale ]
#   -4: Hertzian rectangular contact, curv+half width prescribed,  params = [ mx, my, a1 , bb , scale ]
#   -3: Hertzian elliptical contact, semi-axes prescribed,         params = [ mx, my, aa , bb , scale ]
#   -2: Hertzian elliptical contact, ellipticity prescribed,       params = [ mx, my, a1 , aob, scale ]
#   -1: Hertzian elliptical contact, curvatures prescribed,        params = [ mx, my, a1 , b1 , scale ]
#
# dimensions:  mx, my, aob, scale [-],    a1, b1: [1/length],    aa, bneg, bpos, bb: [length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 1: m=3, cp     - require icp>0, default 1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_sethertzcontact(ire, icp, ipotcn, params):

    if (not isinstance(ire, int)):
        ire = 1
    if (not isinstance(icp, int)):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_sethertzcontact: not available for icp=%d' % icp)

    # convert params to NumPy ndarray
    if (not isinstance(params, np.ndarray)):
        if (isinstance(params, list)):
            params = np.array( params, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_sethertzcontact: invalid params provided.')

    nparam = len(params)

    cntc_dll.cntc_sethertzcontact(c_int(ire), c_int(icp), c_int(ipotcn), c_int(nparam), 
                                                                 params.ctypes.data_as(POINTER(c_double)))

# end function cntc_sethertzcontact

#------------------------------------------------------------------------------------------------------------

