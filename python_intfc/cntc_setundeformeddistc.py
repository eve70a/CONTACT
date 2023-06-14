
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setundeformeddistc(ire, icp, ibase, prmudf)
#
# set the undeformed distance function through a formula or by element-wise specification
#
#  ibase          - type of undeformed distance specification
#  nparam         - number of parameters provided
#  prmudf(nparam) - parameters of undef.dist, depending on method that is used
#
#    1: quadratic function            6    params = [b1, b2, b3, b4, b5, b6]
#    2: circular-x, piecewise-lin-y   5+nn params = [nn, xm, rm, y1, dy1], [b(k), k=1..nn]
#    3: quadratic plus two sines      8    params = [b1, b2, b3, b4, b5, b6, b7, b8]
#    9: elementwise specification     npot params = [h(i), i=1..npot] - undeformed distance per elem. [length]
#                                                   note: positive values == separation between profiles
#
# when ibase=9, prmudf may be of size (my,mx) as well, with nparam=mx*my.
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 1: m=3, cp     - require icp>0, default 1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_setundeformeddistc(ire, icp, ibase, prmudf):

    if (not ire):
        ire = 1
    if (not icp):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_setundeformeddistc: not available for icp=%d' % icp)

    # convert prmudf to NumPy ndarray
    if (not isinstance(prmudf, np.ndarray)):
        if (isinstance(prmudf, list)):
            prmudf = np.array( prmudf, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_setundeformeddistc: invalid prmudf provided.')

    nparam = len(prmudf)

    # TODO: support 2-d size (my,mx), reshape to 1-d array

    cntc_dll.cntc_setundeformeddistc(c_int(ire), c_int(icp), c_int(ibase), c_int(nparam),
                                                                 prmudf.ctypes.data_as(POINTER(c_double)))

# end function cntc_setundeformeddistc

#------------------------------------------------------------------------------------------------------------

