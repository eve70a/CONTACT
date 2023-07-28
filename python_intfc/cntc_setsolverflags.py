
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setsolverflags(ire, icp, imeth, iparam, rparam)
#
# set parameters for the iterative solution algorithms
#  imeth                          - G-digit, which solvers to use
#  nints, nreals                  - number of integer and real parameters provided
#  iparams(nints), rparam(nreals) - depending on method that is used
#
#  0: use default solvers      - iparam = [maxgs, maxin, maxnr, maxout],         rparam = [eps]
#  1: keep parameters          - iparam = [ ],                                   rparam = [ ]
#  2: always use ConvexGS      - iparam = [maxgs, maxin, maxnr, maxout, inislp], rparam = [eps, omegah, omegas, omgslp]
#  3: use SteadyGS if possible - iparam = [maxgs, maxin, maxnr, maxout, inislp], rparam = [eps, omegah, omegas, omgslp]
#  4: use default solvers      - iparam = [maxgs, maxin, maxnr, maxout, inislp], rparam = [eps, omgslp]
#  6: flags sensitivities      - iparam = [mxsens],                              rparam = [epsens]
#
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_setsolverflags(ire, icp, imeth, iparam, rparam):
    # default: W/R contact, all patches

    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp = -1

    # convert iparam and rparam to NumPy ndarrays
    if (not isinstance(iparam, np.ndarray)):
        if (isinstance(iparam, list)):
            iparam = np.array( iparam, dtype=c_int )
        else:
            sys.exit('ERROR in cntc_setsolverflags: invalid iparam provided.')
    if (not isinstance(rparam, np.ndarray)):
        if (isinstance(rparam, list)):
            rparam = np.array( rparam, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_setsolverflags: invalid rparam provided.')

    nints  = len(iparam)
    nreals = len(rparam)

    cntc_dll.cntc_setsolverflags(c_int(ire), c_int(icp), c_int(imeth), 
                                 c_int(nints), iparam.ctypes.data_as(POINTER(c_int)),
                                 c_int(nreals), rparam.ctypes.data_as(POINTER(c_double)))

# end function cntc_setsolverflags

#------------------------------------------------------------------------------------------------------------

