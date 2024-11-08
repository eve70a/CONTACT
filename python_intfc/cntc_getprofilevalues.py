#------------------------------------------------------------------------------------------------------------
# function [ values ] = cntc_getprofilevalues(ire, itask, iparam, rparam)
#
# Get a wheel or rail profile for a wheel-rail contact problem as a table of values
#
#  itask          - select type of outputs:
#                    -1: ierr   <0: error codes, 0=profile loaded ok, 1=not set
#                     0: npnt   number of points used in requested sampling method
#                     1: r/w    get (yr,zr) values for rail or (yw,zw) for wheel profile
#                     2: trk    get (ytr,ztr) values for rail or wheel profile (principal profile)
#                     3: gaug   get left-most point within gauge height, offset, point at gauge height
#                               [ygauge1, zgauge1, yoffs, zoffs, ygauge2, zgauge2]  [length]
#                     4: arc    get arc-length parameter s along profile
#                     5: angl   get surface inclination atan2(dz, dy) [angle]
#                     6: 2d-u   evaluate varprof [xyz]_r or [xyz]_w for cross-section at constant u_i = x_i
#                     7: 2d-v   evaluate varprof [xyz]_r or [xyz]_w for interpolation path at constant v_j
#                     8: 2d-y   evaluate varprof [xyz]_r or [xyz]_w for longitudinal slice at constant y_j
#                     9: spl2d  dump 2D spline data-structure for inspection in Matlab
#  iparam         - integer configuration parameters
#                     1: itype     0 = rail, 1 = wheel profile
#                     2: isampl   -1 = sampling cf. original input data;
#                                  0 = sampling cf. spline representation (default for tasks 0--5);
#                                  1 = sampling cf. spline representation at spacing ds_out (required for
#                                      tasks 6--8)
#                            kchk>=2 = sampling cf. spline representation with integer refinement factor
#  rparam         - real configuration parameters
#                     1: ds_out  step-size ds/du/dv used with sampling method isampl=1, default 1mm
#                     2: c1_out  sample position on variable profiles (task 1, 2: x_out)
#                     3: c2_sta  start of sample range on variable profiles
#                     4: c2_end  end of sample range on variable profiles
#                   task 6: c1 = u, c2 = v; task 7: c1 = v, c2 = u; task 8: c1 = y, c2 = u.
#
#  units: s,x,y,z [mm], th [rad], u,v [-]
#  tasks 1--9: no unit conversion or scaling are applied for profile values
#------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------
# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getprofilevalues(ire, itask, iparam, rparam=None):
    CNTC_err_profil = -32

    if (not isinstance(ire, int)):
        ire    = 1
    if (not isinstance(itask, int)):
        itask  = 0
    if (not iparam):
        itype  = 0       # 0 = rail
        isampl = 0       # 0 = spline sampling
        iparam = [ itype, isampl ]
    if (not rparam):
        ds_out = 1.0     # step size ds
        rparam = [ ds_out ]
    nints  = len(iparam)
    nreals = len(rparam)

    if (itask<-1 or itask>9):
        print('ERROR in cntc_getprofilevalues: itask=%d not available' % itask)
        return -1
    if (not isinstance(iparam, np.ndarray)):
        if (isinstance(iparam, list)):
            iparam = np.array( iparam, dtype=c_int )
        else:
            sys.exit('ERROR in cntc_getprofilevalues: invalid iparam provided.')
    if (not isinstance(rparam, np.ndarray)):
        if (isinstance(rparam, list)):
            rparam = np.array( rparam, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_getprofilevalues: invalid rparam provided.')

    if (itask==-1): # task -1: get error code

        lenarr = 2
        tmp    = np.zeros(lenarr, dtype=c_double)

        cntc_dll.cntc_getprofilevalues_new(c_int(ire), c_int(itask),
                                   c_int(nints),  iparam.ctypes.data_as(POINTER(c_int)),
                                   c_int(nreals), rparam.ctypes.data_as(POINTER(c_double)),
                                   c_int(lenarr), tmp.ctypes.data_as(POINTER(c_double)))

        val    = np.array(tmp, dtype=c_int)

    else:

        # task 0--9: first get the number of points in the profile at requested s sampling positions

        lenarr = 1
        tmp    = np.zeros(lenarr, dtype=c_double)
        if (itask<=5):
            itask_npnt = 0
        else:
            itask_npnt = itask

        cntc_dll.cntc_getprofilevalues_new(c_int(ire), c_int(itask_npnt),
                                c_int(nints),  iparam.ctypes.data_as(POINTER(c_int)),
                                c_int(nreals), rparam.ctypes.data_as(POINTER(c_double)),
                                c_int(lenarr), tmp.ctypes.data_as(POINTER(c_double)))
        npnt = int(tmp[0])

        if (npnt==CNTC_err_profil or npnt<0):
            print('cntc_getprofilevalues: the rail and/or wheel profiles could not be found or processed',
                '(%d).' % npnt)
            return npnt

        # create output-array depending task to be performed

        if (itask==0):

            val = npnt

        elif (itask>=1 and itask<=9):

            lengths = [ 1, npnt*2, npnt*2, 6, npnt, npnt, npnt*3, npnt*3, npnt*3, 1 ]
            lenarr  = lengths[itask]
            val     = np.zeros(lenarr, dtype=c_double)

            cntc_dll.cntc_getprofilevalues_new(c_int(ire), c_int(itask),
                                    c_int(nints),  iparam.ctypes.data_as(POINTER(c_int)),
                                    c_int(nreals), rparam.ctypes.data_as(POINTER(c_double)),
                                    c_int(lenarr), val.ctypes.data_as(POINTER(c_double)))

            if (itask==1 or itask==2):
                val = np.reshape(val, (npnt, 2), order='F')
            elif (itask==6 or itask==7 or itask==8):
                val = np.reshape(val, (npnt, 3), order='F')

        # endif (itask==0 | 1--9)
    # endif (itask==-1)

    return val

# end function cntc_getprofilevalues

#------------------------------------------------------------------------------------------------------------
