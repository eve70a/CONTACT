
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setprofileinputfname(ire, fname, iparam, rparam)
#
# set a wheel or rail profile filename for a wheel-rail contact problem
#
#  fname          - string: name of profile file
#  iparam         - integer configuration parameters
#                     1: itype     0 = rail, 1 = wheel profile, -1 = taken from file extension (default)
#                     2:  -        not used
#                     3: mirrory   0 or -1 = no mirroring (default), 1 = mirror y coordinate values
#                     4: mirrorz   0 = autodetect (default), -1 = no mirroring, 1 = mirror z values
#                     5: errhndl   configuration of error handling. 
#                                   -2 = continue as much as possible, suppress error messages; 
#                                   -1 = suppress warnings; 0: warn and continue (default);
#                                    1 = signal errors and abort
#                     6: ismooth   selection of smoothing method. 0 = original smoothing spline (default),
#                                    1 = weighted PP smoothing spline, 2 = weighted smoothing B-spline (best),
#                                    3 = arc-fit approach
#  rparam         - real configuration parameters
#                     1: sclfac    scaling factor for conversion to [mm], e.g. 1e3 for data given in [m]
#                                  default (sclfac<=0): using the active unit convention
#                     2: smooth    smoothing parameter lambda for non-weighted spline or l_filt for
#                                  weighted splines/arc-fit smoothing
#                     3: maxomit   fraction: signal error if more than maxomit of profile points are
#                                  discarded after cleanup of profile. Default 0.5, use 1 to disable check.
#                     4: zigthrs   angle threshold for zig-zag detection. Default 5/6*pi, >=pi to disable.
#                     5: kinkhigh  angle threshold for kink detection. Default pi/6, >=pi to disable.
#                     6: kinklow   angle threshold for neighbouring points in kink detection. 
#                                  default kinkhigh/5.
#                     7: kinkwid   half-width of window used for kink detection, [len], default 2 mm
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_setprofileinputfname(ire, fname, iparam, rparam):

    if (not isinstance(ire, int)):
        ire    = 1
    if (not iparam):
        iparam = [ -1 ]
    if (not rparam):
        rparam = [ -1., 0. ]

    # convert iparam and rparam to NumPy ndarrays
    if (not isinstance(iparam, np.ndarray)):
        if (isinstance(iparam, list)):
            iparam = np.array( iparam, dtype=c_int )
        else:
            sys.exit('ERROR in cntc_setprofileinputfname: invalid iparam provided.')
    if (not isinstance(rparam, np.ndarray)):
        if (isinstance(rparam, list)):
            rparam = np.array( rparam, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_setprofileinputfname: invalid rparam provided.')

    fname_bytes = fname.encode( 'utf-8' )
    len_fname   = len(fname_bytes)
    nints       = len(iparam)
    nreals      = len(rparam)

    cntc_dll.cntc_setprofileinputfname(c_int(ire), fname_bytes, c_int(len_fname),
                                 c_int(nints), iparam.ctypes.data_as(POINTER(c_int)),
                                 c_int(nreals), rparam.ctypes.data_as(POINTER(c_double)))

# end function cntc_setprofileinputfname

#------------------------------------------------------------------------------------------------------------

