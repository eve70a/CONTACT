
#------------------------------------------------------------------------------------------------------------
# function [ ] = subs_addblock(ire, icp, iblk, isubs, xparam, yparam, zparam)
#
# set the parameters for a block of points for the subsurface stress calculation for a contact problem
#  iblk           - block number; all blocks with this number or higher will be discarded
#  isubs          - type of block specification
#  x/y/zparam     - parameters describing x/y/z-coordinates of the block
#
#  isubs
#    1: xparam = [                 ]  yparam = [                 ]  zparam = [ NZ, ZL, DZ   ]
#    2: xparam = [ IXL, IXINC, IXH ]  yparam = [ IYL, IYINC, IYH ]  zparam = [ NZ, ZL, DZ   ]
#    3: xparam = [ IX(i), i=1:nx   ]  yparam = [ IY(j), j=1:ny   ]  zparam = [ NZ, ZL, DZ   ]
#    5: xparam = [                 ]  yparam = [                 ]  zparam = [ Z(k), k=1:nz ]
#    6: xparam = [ IXL, IXINC, IXH ]  yparam = [ IYL, IYINC, IYH ]  zparam = [ Z(k), k=1:nz ]
#    7: xparam = [ IX(i), i=1:nx   ]  yparam = [ IY(j), j=1:ny   ]  zparam = [ Z(k), k=1:nz ]
#    9: xparam = [ X(i),  i=1:nx   ]  yparam = [ Y(j),  j=1:ny   ]  zparam = [ Z(k), k=1:nz ]
#       ix, iy [-],    x, y, z, zl, dz [length]
#
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 7: m=*, wtd or cp - default icp=-1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def subs_addblock(ire, icp, iblk, isubs, xparam, yparam, zparam):
    # default: fall-back for all contact patches

    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp = -1
    if (not isinstance(iblk, int)):
        iblk = 1

    if (not isinstance(isubs, int) or not isinstance(zparam, int)):
        sys.exit('ERROR in subs_addblock: isubs, zparam are mandatory.')

    if (isubs in [1, 5]):
        xparam = np.array( [ -999 ], dtype=c_double)
        yparam = np.array( [ -999 ], dtype=c_double)
    elif (not xparam or not yparam):
        sys.exit('ERROR in subs_addblock: xparam, yparam are mandatory for isubs=%d.' % isubs)

    # convert xparam/yparam/zparam to NumPy ndarrays
    if (not isinstance(xparam, np.ndarray)):
        xparam = np.array( xparam, dtype=c_double )
    if (not isinstance(yparam, np.ndarray)):
        yparam = np.array( yparam, dtype=c_double )
    if (not isinstance(zparam, np.ndarray)):
        zparam = np.array( zparam, dtype=c_double )

    npx = len(xparam)
    npy = len(yparam)
    npz = len(zparam)

    cntc_dll.subs_addblock(c_int(ire), c_int(icp), c_int(iblk), c_int(isubs),
                           c_int(npx), c_int(npy), c_int(npz), xparam.ctypes.data_as(POINTER(c_double)),
                           yparam.ctypes.data_as(POINTER(c_double)), zparam.ctypes.data_as(POINTER(c_double)))

# end function subs_addblock

#------------------------------------------------------------------------------------------------------------

