#------------------------------------------------------------------------------------------------------------
# function [ values ] = cntc_getprofilevalues(ire, itask, iparam)
#
# Get a wheel or rail profile for a wheel-rail contact problem as a table of values
#
#  itask          - select type of outputs:
#                     0: npnt   number of points used in profile
#                     1: r/w    get (yr,zr) value for rail or (yw,zw) for wheel profile
#                     2: trk    get (ytr,ztr) values for rail or wheel profile (principal profile)
#                     3: gaug   get left-most point within gauge height, offsets, point at gauge height
#                               [ygauge1, zgauge1, yoffs, zoffs, ygauge2, zgauge2]  [length]
#                     4: arc    get arc-length parameter s along profile
#                     5: angl   get surface inclination atan2(dz, dy) [length]
#  iparam         - integer configuration parameters
#                     1: itype     0 = rail, 1 = wheel profile
#                     2: iside     0 = left, 1 = right side
#  tasks 1,2,4: no unit conversion or scaling are applied for profile values
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getprofilevalues(ire, itask, iparam):
    CNTC_err_profil = -32

    if (not ire):
        ire = 1
    if (not itask):
        itask = 0
    if (not iparam):
        itype = 0       # 0 = rail
        iside = 1       # 1 = right side
        iparam = [ itype, iside ]
    nints = len(iparam)

    if (itask<0 or itask>5):
        print('ERROR in cntc_getprofilevalues: itask=%d not available' % itask)
        return -1
    if (not isinstance(iparam, np.ndarray)):
        if (isinstance(iparam, list)):
            iparam = np.array( iparam, dtype=c_int )
        else:
            sys.exit('ERROR in cntc_setprofileinputvalues: invalid iparam provided.')

    # always get the number of points in the profile (task=0)

    lenarr = 1
    tmp    = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getprofilevalues(c_int(ire), c_int(0),
                               c_int(nints),  iparam.ctypes.data_as(POINTER(c_int)),
                               c_int(lenarr), tmp.ctypes.data_as(POINTER(c_double)))
    npnt = int(tmp[0])

    if (npnt==CNTC_err_profil or npnt<0):
        print('cntc_getprofilevalues: the rail and/or wheel profiles could not be found or processed',
              '(%d).' % npnt)
        return npnt

    # create output-array depending task to be performed

    if (itask==0):

        val = npnt

    elif (itask>=1 and itask<=5):

        lengths = [ 1, npnt*2, npnt*2, 6, npnt, npnt ]
        lenarr  = lengths[itask]
        val     = np.zeros(lenarr, dtype=c_double)

        cntc_dll.cntc_getprofilevalues(c_int(ire), c_int(itask),
                                c_int(nints),  iparam.ctypes.data_as(POINTER(c_int)),
                                c_int(lenarr), val.ctypes.data_as(POINTER(c_double)))

        if (itask==1 or itask==2):
            val = np.reshape(val, (npnt, 2), order='F')

    return val

# end function cntc_getprofilevalues

#------------------------------------------------------------------------------------------------------------

